//! A minimal HTTP/1.1 client over `std::net`, for downloading public data files (REBASE serves
//! plain `http://`). GET only, no TLS, no cookies, no personal data: the request carries just
//! the host, path and the generic `SnappieSnap-dev` user agent (AGENTS.md). Follows up to a few `http://`
//! redirects, decodes chunked bodies and caps size and time.

use std::io::{self, Read, Write};
use std::net::{TcpStream, ToSocketAddrs};
use std::time::{Duration, Instant};

/// Largest body we accept.
pub const MAX_DOWNLOAD: usize = 64 << 20;
const MAX_REDIRECTS: usize = 3;
const MAX_HEADER: usize = 64 << 10;
const TIMEOUT: Duration = Duration::from_secs(30);
const DEADLINE: Duration = Duration::from_secs(300);

fn invalid(msg: impl Into<String>) -> io::Error {
    io::Error::new(io::ErrorKind::InvalidData, msg.into())
}

/// `http://host[:port]/path` → (host, port, path).
fn parse_url(url: &str) -> io::Result<(String, u16, String)> {
    let rest = url.trim().strip_prefix("http://").ok_or_else(|| {
        io::Error::new(
            io::ErrorKind::Unsupported,
            format!("{url}: only http:// URLs can be downloaded"),
        )
    })?;
    let (authority, path) = match rest.find('/') {
        Some(i) => rest.split_at(i),
        None => (rest, "/"),
    };
    if authority.is_empty() || authority.contains(['@', ' ']) || path.contains([' ', '\r', '\n']) {
        return Err(invalid(format!("{url}: not a valid URL")));
    }
    let (host, port) = match authority.rsplit_once(':') {
        Some((h, p)) => (
            h,
            p.parse().map_err(|_| invalid(format!("{url}: bad port")))?,
        ),
        None => (authority, 80),
    };
    Ok((host.to_string(), port, path.to_string()))
}

/// GET `url` and return the body of a 200 response.
pub fn get(url: &str) -> io::Result<Vec<u8>> {
    let started = Instant::now();
    let mut url = url.to_string();
    for _ in 0..=MAX_REDIRECTS {
        let (host, port, path) = parse_url(&url)?;
        match request(&host, port, &path, started)? {
            Response::Body(b) => return Ok(b),
            Response::Redirect(to) => {
                url = if to.starts_with('/') {
                    format!("http://{host}:{port}{to}")
                } else {
                    to
                };
            }
        }
    }
    Err(invalid("too many redirects"))
}

enum Response {
    Body(Vec<u8>),
    Redirect(String),
}

fn request(host: &str, port: u16, path: &str, started: Instant) -> io::Result<Response> {
    let addrs: Vec<_> = (host, port).to_socket_addrs()?.collect();
    let mut stream = None;
    let mut last = io::Error::new(io::ErrorKind::NotFound, format!("{host}: no address"));
    for addr in addrs {
        match TcpStream::connect_timeout(&addr, TIMEOUT) {
            Ok(s) => {
                stream = Some(s);
                break;
            }
            Err(e) => last = e,
        }
    }
    let mut stream = stream.ok_or(last)?;
    stream.set_read_timeout(Some(TIMEOUT))?;
    stream.set_write_timeout(Some(TIMEOUT))?;
    let host_header = if port == 80 {
        host.to_string()
    } else {
        format!("{host}:{port}")
    };
    write!(
        stream,
        "GET {path} HTTP/1.1\r\nHost: {host_header}\r\nUser-Agent: SnappieSnap-dev\r\nAccept: */*\r\nAccept-Encoding: identity\r\nConnection: close\r\n\r\n"
    )?;
    let mut raw = Vec::new();
    let mut buf = [0u8; 64 << 10];
    loop {
        if started.elapsed() > DEADLINE {
            return Err(io::Error::new(
                io::ErrorKind::TimedOut,
                "download took too long",
            ));
        }
        let n = stream.read(&mut buf)?;
        if n == 0 {
            break;
        }
        raw.extend_from_slice(buf.get(..n).unwrap_or_default());
        if raw.len() > MAX_DOWNLOAD + MAX_HEADER {
            return Err(invalid(format!(
                "download larger than {} MiB",
                MAX_DOWNLOAD >> 20
            )));
        }
    }
    parse_response(&raw)
}

fn parse_response(raw: &[u8]) -> io::Result<Response> {
    let split = raw
        .windows(4)
        .position(|w| w == b"\r\n\r\n")
        .ok_or_else(|| invalid("incomplete HTTP response"))?;
    let (head, body) = raw.split_at(split);
    let body = body.get(4..).unwrap_or_default();
    let head = String::from_utf8_lossy(head);
    let mut lines = head.split("\r\n");
    let status_line = lines.next().unwrap_or_default();
    let status: u16 = status_line
        .split_whitespace()
        .nth(1)
        .and_then(|s| s.parse().ok())
        .ok_or_else(|| invalid(format!("bad HTTP status line: {status_line}")))?;
    let mut length = None;
    let mut chunked = false;
    let mut location = None;
    for line in lines {
        let Some((name, value)) = line.split_once(':') else {
            continue;
        };
        let value = value.trim();
        match name.trim().to_ascii_lowercase().as_str() {
            "content-length" => length = value.parse::<usize>().ok(),
            "transfer-encoding" => chunked = value.to_ascii_lowercase().contains("chunked"),
            "location" => location = Some(value.to_string()),
            _ => {}
        }
    }
    match status {
        200 => {}
        301 | 302 | 303 | 307 | 308 => {
            return location
                .map(Response::Redirect)
                .ok_or_else(|| invalid(format!("HTTP {status} without a Location")));
        }
        _ => return Err(io::Error::other(format!("server answered: {status_line}"))),
    }
    let body = if chunked {
        dechunk(body)?
    } else {
        match length {
            Some(n) => body
                .get(..n)
                .ok_or_else(|| invalid("connection closed before the download finished"))?
                .to_vec(),
            None => body.to_vec(),
        }
    };
    Ok(Response::Body(body))
}

fn dechunk(mut body: &[u8]) -> io::Result<Vec<u8>> {
    let mut out = Vec::new();
    loop {
        let eol = body
            .windows(2)
            .position(|w| w == b"\r\n")
            .ok_or_else(|| invalid("bad chunked encoding"))?;
        let size_line = String::from_utf8_lossy(body.get(..eol).unwrap_or_default()).into_owned();
        let size_hex = size_line.split(';').next().unwrap_or_default().trim();
        let size = usize::from_str_radix(size_hex, 16).map_err(|_| invalid("bad chunk size"))?;
        let start = eol + 2;
        if size == 0 {
            return Ok(out);
        }
        let end = start
            .checked_add(size)
            .filter(|&e| e <= body.len())
            .ok_or_else(|| invalid("truncated chunk"))?;
        out.extend_from_slice(body.get(start..end).unwrap_or_default());
        if out.len() > MAX_DOWNLOAD {
            return Err(invalid("download too large"));
        }
        body = body.get(end + 2..).unwrap_or_default();
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn urls() {
        assert_eq!(
            parse_url("http://rebase.neb.com/rebase/link_withrefm").unwrap(),
            ("rebase.neb.com".into(), 80, "/rebase/link_withrefm".into())
        );
        assert_eq!(
            parse_url("http://h:8080").unwrap(),
            ("h".into(), 8080, "/".into())
        );
        for bad in [
            "https://h/x",
            "ftp://h",
            "http://",
            "http://u@h/",
            "http://h:x/",
            "http://h/a b",
        ] {
            assert!(parse_url(bad).is_err(), "{bad}");
        }
    }

    #[test]
    fn responses() {
        let ok = parse_response(b"HTTP/1.1 200 OK\r\nContent-Length: 3\r\n\r\nabcdef").unwrap();
        assert!(matches!(ok, Response::Body(b) if b == b"abc"));
        let chunked =
            parse_response(b"HTTP/1.1 200 OK\r\nTransfer-Encoding: chunked\r\n\r\n3\r\nabc\r\n2;x\r\nde\r\n0\r\n\r\n")
                .unwrap();
        assert!(matches!(chunked, Response::Body(b) if b == b"abcde"));
        let moved = parse_response(b"HTTP/1.1 301 Moved\r\nLocation: http://x/y\r\n\r\n").unwrap();
        assert!(matches!(moved, Response::Redirect(u) if u == "http://x/y"));
        for bad in [
            &b""[..],
            b"HTTP/1.1 404 Not Found\r\n\r\n",
            b"garbage\r\n\r\n",
            b"HTTP/1.1 200 OK\r\nContent-Length: 10\r\n\r\nabc",
            b"HTTP/1.1 200 OK\r\nTransfer-Encoding: chunked\r\n\r\nffffffffffffffff\r\nabc",
            b"HTTP/1.1 200 OK\r\nTransfer-Encoding: chunked\r\n\r\nzz\r\n",
            b"HTTP/1.1 302 Found\r\n\r\n",
        ] {
            assert!(
                parse_response(bad).is_err(),
                "{}",
                String::from_utf8_lossy(bad)
            );
        }
    }
}
