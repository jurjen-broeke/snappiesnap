//! A small XML tokenizer for the simple, attribute-heavy XML inside `.dna` packets.
//!
//! It reports each token with its byte span so callers can keep raw text, never recurses, and
//! returns an error (never panics) on malformed input.

use std::ops::Range;

#[derive(Debug, Clone, PartialEq, Eq, thiserror::Error)]
#[error("malformed XML at byte {offset}: {reason}")]
pub struct XmlError {
    pub offset: usize,
    pub reason: &'static str,
}

/// An attribute with its value exactly as written (still escaped).
#[derive(Debug, Clone, PartialEq, Eq)]
pub struct RawAttr<'a> {
    pub name: &'a str,
    pub raw: &'a str,
}

#[derive(Debug, Clone, PartialEq, Eq)]
pub enum Token<'a> {
    Start {
        name: &'a str,
        attrs: Vec<RawAttr<'a>>,
        self_closing: bool,
    },
    End {
        name: &'a str,
    },
    Text(&'a str),
    /// `<?...?>`, `<!--...-->` or `<!...>`.
    Other(&'a str),
}

pub struct Tokenizer<'a> {
    src: &'a str,
    pos: usize,
}

impl<'a> Tokenizer<'a> {
    pub fn new(src: &'a str) -> Self {
        Tokenizer { src, pos: 0 }
    }

    fn err(&self, offset: usize, reason: &'static str) -> XmlError {
        let _ = self;
        XmlError { offset, reason }
    }

    fn rest(&self) -> &'a str {
        self.src.get(self.pos..).unwrap_or("")
    }

    /// The next token and its span, `None` at the end.
    pub fn next_token(&mut self) -> Option<Result<(Token<'a>, Range<usize>), XmlError>> {
        let start = self.pos;
        let rest = self.rest();
        if rest.is_empty() {
            return None;
        }
        let result = if !rest.starts_with('<') {
            let len = rest.find('<').unwrap_or(rest.len());
            self.pos = start.saturating_add(len);
            Ok(Token::Text(rest.get(..len).unwrap_or("")))
        } else if rest.starts_with("<?") {
            self.until(start, "?>")
        } else if rest.starts_with("<!--") {
            self.until(start, "-->")
        } else if rest.starts_with("<!") {
            self.until(start, ">")
        } else if rest.starts_with("</") {
            self.end_tag(start)
        } else {
            self.start_tag(start)
        };
        Some(result.map(|t| (t, start..self.pos)))
    }

    fn until(&mut self, start: usize, terminator: &str) -> Result<Token<'a>, XmlError> {
        let rest = self.rest();
        let end = rest
            .find(terminator)
            .ok_or_else(|| self.err(start, "unterminated declaration or comment"))?
            .saturating_add(terminator.len());
        self.pos = start.saturating_add(end);
        Ok(Token::Other(rest.get(..end).unwrap_or("")))
    }

    fn end_tag(&mut self, start: usize) -> Result<Token<'a>, XmlError> {
        let rest = self.rest();
        let close = rest
            .find('>')
            .ok_or_else(|| self.err(start, "unterminated end tag"))?;
        let name = rest.get(2..close).unwrap_or("").trim();
        if !valid_name(name) {
            return Err(self.err(start, "invalid end tag name"));
        }
        self.pos = start.saturating_add(close).saturating_add(1);
        Ok(Token::End { name })
    }

    fn start_tag(&mut self, start: usize) -> Result<Token<'a>, XmlError> {
        let src = self.src;
        let bytes = src.as_bytes();
        let mut i = start.saturating_add(1);
        let name_start = i;
        while bytes
            .get(i)
            .is_some_and(|b| !b.is_ascii_whitespace() && *b != b'>' && *b != b'/')
        {
            i = i.saturating_add(1);
        }
        let name = src.get(name_start..i).unwrap_or("");
        if !valid_name(name) {
            return Err(self.err(start, "invalid tag name"));
        }
        let mut attrs = Vec::new();
        loop {
            while bytes.get(i).is_some_and(u8::is_ascii_whitespace) {
                i = i.saturating_add(1);
            }
            match bytes.get(i) {
                None => return Err(self.err(start, "unterminated start tag")),
                Some(b'>') => {
                    self.pos = i.saturating_add(1);
                    return Ok(Token::Start {
                        name,
                        attrs,
                        self_closing: false,
                    });
                }
                Some(b'/') => {
                    if bytes.get(i.saturating_add(1)) != Some(&b'>') {
                        return Err(self.err(i, "expected '>' after '/'"));
                    }
                    self.pos = i.saturating_add(2);
                    return Ok(Token::Start {
                        name,
                        attrs,
                        self_closing: true,
                    });
                }
                Some(_) => {}
            }
            let attr_start = i;
            while bytes
                .get(i)
                .is_some_and(|b| *b != b'=' && !b.is_ascii_whitespace() && *b != b'>')
            {
                i = i.saturating_add(1);
            }
            let attr_name = src.get(attr_start..i).unwrap_or("");
            while bytes.get(i).is_some_and(u8::is_ascii_whitespace) {
                i = i.saturating_add(1);
            }
            if bytes.get(i) != Some(&b'=') || !valid_name(attr_name) {
                return Err(self.err(attr_start, "expected name=\"value\""));
            }
            i = i.saturating_add(1);
            while bytes.get(i).is_some_and(u8::is_ascii_whitespace) {
                i = i.saturating_add(1);
            }
            let quote = match bytes.get(i) {
                Some(q @ (b'"' | b'\'')) => *q,
                _ => return Err(self.err(i, "attribute value must be quoted")),
            };
            let value_start = i.saturating_add(1);
            let value_len = src
                .get(value_start..)
                .and_then(|r| r.as_bytes().iter().position(|b| *b == quote))
                .ok_or_else(|| self.err(i, "unterminated attribute value"))?;
            let value_end = value_start.saturating_add(value_len);
            attrs.push(RawAttr {
                name: attr_name,
                raw: src.get(value_start..value_end).unwrap_or(""),
            });
            i = value_end.saturating_add(1);
        }
    }
}

fn valid_name(name: &str) -> bool {
    !name.is_empty()
        && name
            .bytes()
            .all(|b| b.is_ascii_alphanumeric() || matches!(b, b'_' | b'-' | b'.' | b':'))
}

/// Replace the five predefined entities and numeric character references. Unknown or
/// malformed references are kept literally.
pub fn decode(raw: &str) -> String {
    if !raw.contains('&') {
        return raw.to_string();
    }
    let mut out = String::with_capacity(raw.len());
    let mut rest = raw;
    while let Some(amp) = rest.find('&') {
        out.push_str(rest.get(..amp).unwrap_or(""));
        let after = rest.get(amp..).unwrap_or("");
        let entity_end = after.find(';').filter(|&e| e <= 12);
        let replacement = entity_end.and_then(|end| {
            let name = after.get(1..end)?;
            let c = match name {
                "lt" => '<',
                "gt" => '>',
                "amp" => '&',
                "quot" => '"',
                "apos" => '\'',
                _ => {
                    let code = if let Some(hex) = name.strip_prefix("#x") {
                        u32::from_str_radix(hex, 16).ok()?
                    } else {
                        name.strip_prefix('#')?.parse::<u32>().ok()?
                    };
                    char::from_u32(code)?
                }
            };
            Some((c, end))
        });
        match replacement {
            Some((c, end)) => {
                out.push(c);
                rest = after.get(end.saturating_add(1)..).unwrap_or("");
            }
            None => {
                out.push('&');
                rest = after.get(1..).unwrap_or("");
            }
        }
    }
    out.push_str(rest);
    out
}

/// Escape an attribute value the way SnapGene writes feature XML: `& < > "` become entities,
/// line breaks and tabs become character references (`&#10;`, as observed; raw ones would be
/// normalised to spaces by XML readers); everything else (including `'` and non-ASCII) is
/// written as is.
pub fn encode_attr(value: &str) -> String {
    let mut out = String::with_capacity(value.len());
    for c in value.chars() {
        match c {
            '&' => out.push_str("&amp;"),
            '<' => out.push_str("&lt;"),
            '>' => out.push_str("&gt;"),
            '"' => out.push_str("&quot;"),
            '\n' => out.push_str("&#10;"),
            '\r' => out.push_str("&#13;"),
            '\t' => out.push_str("&#9;"),
            other => out.push(other),
        }
    }
    out
}

#[cfg(test)]
mod tests {
    use super::*;

    fn tokens(src: &str) -> Result<Vec<Token<'_>>, XmlError> {
        let mut t = Tokenizer::new(src);
        let mut out = Vec::new();
        while let Some(tok) = t.next_token() {
            out.push(tok?.0);
        }
        Ok(out)
    }

    #[test]
    fn tokenizes_snapgene_style_xml() {
        let toks =
            tokens(r#"<?xml version="1.0"?><A x="1" y='2'><B z="&lt;i&gt;"/>text</A>"#).unwrap();
        assert_eq!(toks.len(), 5);
        assert_eq!(toks[0], Token::Other(r#"<?xml version="1.0"?>"#));
        match &toks[2] {
            Token::Start {
                name,
                attrs,
                self_closing,
            } => {
                assert_eq!(*name, "B");
                assert!(*self_closing);
                assert_eq!(attrs[0].raw, "&lt;i&gt;");
            }
            other => panic!("{other:?}"),
        }
        assert_eq!(toks[3], Token::Text("text"));
    }

    #[test]
    fn malformed_input_is_an_error() {
        for bad in [
            "<A",
            "<A x>",
            "<A x=1>",
            "<A x=\"1>",
            "</A",
            "<?x",
            "<!--",
            "<A/x>",
            "< A>",
            "<>",
        ] {
            assert!(tokens(bad).is_err(), "{bad:?}");
        }
    }

    #[test]
    fn entities() {
        assert_eq!(
            decode("a &lt;b&gt; &amp; &quot;c&quot; &apos;"),
            "a <b> & \"c\" '"
        );
        assert_eq!(
            decode("&#169; &#xAE; &bogus; & &#xFFFFFFFF;"),
            "© ® &bogus; & &#xFFFFFFFF;"
        );
        let s = "x < y & \"z\" 'q' ®";
        assert_eq!(decode(&encode_attr(s)), s);
        assert_eq!(encode_attr("aph(3')"), "aph(3')");
        assert_eq!(encode_attr("line\nbreak\t"), "line&#10;break&#9;");
        assert_eq!(decode("line&#10;break"), "line\nbreak");
    }
}
