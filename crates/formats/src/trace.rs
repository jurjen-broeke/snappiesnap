//! A sequencing trace (chromatogram) as the trace formats ([`crate::abif`], [`crate::scf`])
//! describe it: four signal channels, the called bases, where each base's peak is in the
//! signal, and a quality per base.

/// The channels in the order [`Chromatogram::channels`] holds them.
pub const CHANNEL_BASES: [u8; 4] = *b"ACGT";

/// A chromatogram read from a trace file.
#[derive(Debug, Clone, Default, PartialEq, Eq)]
pub struct Chromatogram {
    /// The sample name recorded by the instrument, if any.
    pub sample_name: Option<String>,
    /// The called bases, as stored (IUPAC letters; case kept).
    pub bases: Vec<u8>,
    /// Per base: the index of its peak in the channels.
    pub peaks: Vec<u32>,
    /// Per base: its quality on the Phred scale (`-10 log10 P(error)`), if the file has them.
    pub quality: Option<Vec<u8>>,
    /// The signal, one value per sample point, for A, C, G and T (see [`CHANNEL_BASES`]). All
    /// four have the same length.
    pub channels: [Vec<u16>; 4],
}

impl Chromatogram {
    /// The number of sample points per channel.
    pub fn samples(&self) -> usize {
        self.channels.first().map_or(0, Vec::len)
    }

    /// Check that the per-base and per-sample arrays agree in length.
    pub fn check(&self) -> Result<(), String> {
        let n = self.bases.len();
        if self.peaks.len() != n {
            return Err(format!("{n} bases but {} peak positions", self.peaks.len()));
        }
        if let Some(q) = &self.quality
            && q.len() != n
        {
            return Err(format!("{n} bases but {} quality values", q.len()));
        }
        let samples = self.samples();
        if self.channels.iter().any(|c| c.len() != samples) {
            return Err("the four signal channels differ in length".into());
        }
        Ok(())
    }
}
