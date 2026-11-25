//! The genetic codes, generated from NCBI's `gc.prt` (version 4.6,
//! https://ftp.ncbi.nih.gov/entrez/misc/data/gc.prt), a work of the US government and in the
//! public domain. Codons are in NCBI's order: first, second and third base each run T, C, A, G.
//! In `starts`, `M` marks an initiation codon and `*` a stop.
//!
//! Generated; edit by regenerating from the source file, not by hand.

/// One genetic code.
#[derive(Debug, PartialEq, Eq)]
pub struct GeneticCode {
    /// NCBI's translation table number.
    pub id: u8,
    pub name: &'static str,
    /// Amino acid (one letter, `*` for stop) of each of the 64 codons.
    pub amino_acids: &'static [u8; 64],
    /// Start-codon flags of each codon (`M` = can initiate).
    pub starts: &'static [u8; 64],
}

/// Every NCBI genetic code, by table number.
pub const CODES: &[GeneticCode] = &[
    GeneticCode {
        id: 1,
        name: "Standard",
        amino_acids: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**--*----M---------------M----------------------------",
    },
    GeneticCode {
        id: 2,
        name: "Vertebrate Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSS**VVVVAAAADDEEGGGG",
        starts: b"----------**--------------------MMMM----------**---M------------",
    },
    GeneticCode {
        id: 3,
        name: "Yeast Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWTTTTPPPPHHQQRRRRIIMMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**----------------------MM---------------M------------",
    },
    GeneticCode {
        id: 4,
        name: "Mold Mitochondrial; Protozoan Mitochondrial; Coelenterate Mitochondrial; Mycoplasma; Spiroplasma",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--MM------**-------M------------MMMM---------------M------------",
    },
    GeneticCode {
        id: 5,
        name: "Invertebrate Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSSSVVVVAAAADDEEGGGG",
        starts: b"---M------**--------------------MMMM---------------M------------",
    },
    GeneticCode {
        id: 6,
        name: "Ciliate Nuclear; Dasycladacean Nuclear; Hexamita Nuclear",
        amino_acids: b"FFLLSSSSYYQQCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    GeneticCode {
        id: 9,
        name: "Echinoderm Mitochondrial; Flatworm Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M---------------M------------",
    },
    GeneticCode {
        id: 10,
        name: "Euplotid Nuclear",
        amino_acids: b"FFLLSSSSYY**CCCWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M----------------------------",
    },
    GeneticCode {
        id: 11,
        name: "Bacterial, Archaeal and Plant Plastid",
        amino_acids: b"FFLLSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**--*----M------------MMMM---------------M------------",
    },
    GeneticCode {
        id: 12,
        name: "Alternative Yeast Nuclear",
        amino_acids: b"FFLLSSSSYY**CC*WLLLSPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*----M---------------M----------------------------",
    },
    GeneticCode {
        id: 13,
        name: "Ascidian Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNKKSSGGVVVVAAAADDEEGGGG",
        starts: b"---M------**----------------------MM---------------M------------",
    },
    GeneticCode {
        id: 14,
        name: "Alternative Flatworm Mitochondrial",
        amino_acids: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"-----------*-----------------------M----------------------------",
    },
    GeneticCode {
        id: 15,
        name: "Blepharisma Macronuclear",
        amino_acids: b"FFLLSSSSYY*QCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------*---*--------------------M----------------------------",
    },
    GeneticCode {
        id: 16,
        name: "Chlorophycean Mitochondrial",
        amino_acids: b"FFLLSSSSYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------*---*--------------------M----------------------------",
    },
    GeneticCode {
        id: 21,
        name: "Trematode Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIMMTTTTNNNKSSSSVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M---------------M------------",
    },
    GeneticCode {
        id: 22,
        name: "Scenedesmus obliquus Mitochondrial",
        amino_acids: b"FFLLSS*SYY*LCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"------*---*---*--------------------M----------------------------",
    },
    GeneticCode {
        id: 23,
        name: "Thraustochytrium Mitochondrial",
        amino_acids: b"FF*LSSSSYY**CC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--*-------**--*-----------------M--M---------------M------------",
    },
    GeneticCode {
        id: 24,
        name: "Rhabdopleuridae Mitochondrial",
        amino_acids: b"FFLLSSSSYY**CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts: b"---M------**-------M---------------M---------------M------------",
    },
    GeneticCode {
        id: 25,
        name: "Candidate Division SR1 and Gracilibacteria",
        amino_acids: b"FFLLSSSSYY**CCGWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------**-----------------------M---------------M------------",
    },
    GeneticCode {
        id: 26,
        name: "Pachysolen tannophilus Nuclear",
        amino_acids: b"FFLLSSSSYY**CC*WLLLAPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*----M---------------M----------------------------",
    },
    GeneticCode {
        id: 27,
        name: "Karyorelict Nuclear",
        amino_acids: b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    GeneticCode {
        id: 28,
        name: "Condylostoma Nuclear",
        amino_acids: b"FFLLSSSSYYQQCCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**--*--------------------M----------------------------",
    },
    GeneticCode {
        id: 29,
        name: "Mesodinium Nuclear",
        amino_acids: b"FFLLSSSSYYYYCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    GeneticCode {
        id: 30,
        name: "Peritrich Nuclear",
        amino_acids: b"FFLLSSSSYYEECC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"--------------*--------------------M----------------------------",
    },
    GeneticCode {
        id: 31,
        name: "Blastocrithidia Nuclear",
        amino_acids: b"FFLLSSSSYYEECCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"----------**-----------------------M----------------------------",
    },
    GeneticCode {
        id: 32,
        name: "Balanophoraceae Plastid",
        amino_acids: b"FFLLSSSSYY*WCC*WLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSRRVVVVAAAADDEEGGGG",
        starts: b"---M------*---*----M------------MMMM---------------M------------",
    },
    GeneticCode {
        id: 33,
        name: "Cephalodiscidae Mitochondrial",
        amino_acids: b"FFLLSSSSYYY*CCWWLLLLPPPPHHQQRRRRIIIMTTTTNNKKSSSKVVVVAAAADDEEGGGG",
        starts: b"---M-------*-------M---------------M---------------M------------",
    },
];
