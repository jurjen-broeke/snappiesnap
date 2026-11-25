use std::collections::VecDeque;

use serde::Serialize;

/// Entries kept in memory; older ones are dropped.
const CAPACITY: usize = 1000;

#[derive(Debug, Clone, Copy, PartialEq, Eq, Serialize)]
#[serde(rename_all = "lowercase")]
pub enum Level {
    Info,
    Warn,
    Error,
}

#[derive(Debug, Clone, PartialEq, Eq, Serialize)]
pub struct LogEntry {
    /// Increases by one per entry, so readers can ask for "everything after n".
    pub seq: u64,
    pub level: Level,
    /// The command or subsystem that logged it.
    pub source: String,
    pub message: String,
}

/// A bounded in-memory log of recent problems (exposed as `log.tail`).
#[derive(Debug, Default, Clone)]
pub(crate) struct Log {
    entries: VecDeque<LogEntry>,
    next_seq: u64,
}

impl Log {
    pub fn push(&mut self, level: Level, source: &str, message: impl Into<String>) {
        if self.entries.len() >= CAPACITY {
            self.entries.pop_front();
        }
        self.entries.push_back(LogEntry {
            seq: self.next_seq,
            level,
            source: source.to_string(),
            message: message.into(),
        });
        self.next_seq = self.next_seq.saturating_add(1);
    }

    /// The last `n` entries, oldest first.
    pub fn tail(&self, n: usize) -> Vec<LogEntry> {
        let skip = self.entries.len().saturating_sub(n);
        self.entries.iter().skip(skip).cloned().collect()
    }
}

#[cfg(test)]
mod tests {
    use super::*;

    #[test]
    fn keeps_only_the_most_recent_entries() {
        let mut log = Log::default();
        for i in 0..CAPACITY + 5 {
            log.push(Level::Info, "t", format!("{i}"));
        }
        let all = log.tail(usize::MAX);
        assert_eq!(all.len(), CAPACITY);
        assert_eq!(all[0].seq, 5);
        assert_eq!(log.tail(2).len(), 2);
        assert!(log.tail(0).is_empty());
    }
}
