# Built-in features: sources

SnappieSnap's built-in features (`insdc_features.json`) are cut from public INSDC/RefSeq
records by `cargo xtask features`, following `manifest.json`. Each entry is the sequence of
one feature at the location given by **that record's own feature table**. The generator checks
that the record has a feature with that key and location whose qualifiers name the expected
element, and that every CDS translates to the record's `/translation`.

## Which records may be used

Only records from genome projects, academic laboratories and NCBI's curated RefSeq.
Records submitted by companies (vector suppliers and others) are **never** used, even when
they are public. Third-party feature databases (including those of other sequence editors)
are never used.

## Terms

NCBI places no restrictions on the use or distribution of GenBank data
(https://www.ncbi.nlm.nih.gov/home/about/policies/). It notes that submitters may claim
rights in their records. This is why only genome, academic and RefSeq records are used.
The names, types and notes in the manifest are our own wording.

## Records

| Accession | Elements | Reference (from the record) |
|---|---|---|
| J01749.1 | AmpR, TcR, rop | Sutcliffe, J.G., Proc. Natl. Acad. Sci. USA 75:3737 (1978): pBR322 |
| J01636.1 | lacI, CAP binding site, lac operator, lacZ | Gilbert, W. & Maxam, A., Proc. Natl. Acad. Sci. USA 70:3581 (1973) and later papers: E. coli lac operon |
| NC_001669.1 | SV40 ori, SV40 72 bp enhancer repeat, SV40 21 bp repeat | NCBI RefSeq, Simian virus 40 genome |
| V00618.1 | NeoR/KanR | Beck, E. et al., Gene 19:327 (1982): Tn5 |
| V00359.1 | KanR (Tn903) | Oka, A., Sugisaki, H. & Takanami, M., J. Mol. Biol. 147:217 (1981) |
| V00622.1 | CmR | Alton, N.K. & Vapnek, D., Nature 282:864 (1979): Tn9 |
| K01193.1 | HygR | Gritz, L. & Davies, J., Gene 25:179 (1983) |
| M25346.1 | PuroR | Lacalle, R.A. et al., Gene 79:375 (1989) |
| X03453.1 | loxP | Sternberg, N. et al., J. Mol. Biol. 187:197 (1986): phage P1 |
| J01347.1 | FLP binding site | Hindley, J. & Phear, G.A., Nucleic Acids Res. 7:361 (1979) and later papers: yeast 2 micron plasmid |
| AY678264.1 | mCherry | Shaner, N.C. et al., Nat. Biotechnol. 22:1567 (2004) |
| M62653.1 | GFP (wild type) | Prasher, D.C. et al., Gene 111:229 (1992) |
| K03104.1 | CMV enhancer | Boshart, M. et al., Cell 41:521 (1985) |
| J01830.1 | TetR | Bertrand, K.P. et al. and Hillen, W. et al. (1982-1984): transposon Tn10 |
| J01641.1 | AraC | Smith, B.R. & Schleif, R., J. Biol. Chem. 253:6931 (1978) |
| NC_001416.1 | lambda cI | NCBI RefSeq, Enterobacteria phage lambda genome |
| NC_003287.2 | M13 gIII | NCBI RefSeq, Enterobacteria phage M13 genome |
| X73149.1 | NatR | Krugel, H. et al., Gene 127:127 (1993) |
| D83710.1 | BSD | Kimura, M. et al. (1994), Aspergillus terreus |
| X02730.1 | sacB | Fouet, A. et al. (1984), Bacillus subtilis |
| NC_002483.1 | ccdA, ccdB | NCBI RefSeq, E. coli K-12 F plasmid |
| V01146.1 | T7 promoter | Dunn, J.J. & Studier, F.W., J. Mol. Biol. 166:477 (1983), PMID 6864790; boundaries from the paper (see below) |
| NC_001137.3, NC_001135.5, NC_001147.6, NC_001136.10 | URA3, LEU2, HIS3, TRP1 | NCBI RefSeq, Saccharomyces cerevisiae S288C (from SGD) |
| LC375816.1 | attB1, attB2 | Sugano, S.S. et al., PLoS ONE 13:e0205117 (2018): binary vector pMpGE014 (Ritsumeikan University) |
| LC090754.1 | attL1, attL2 | Sugano, S.S. et al., PLoS ONE 13:e0205117 (2018): Gateway entry vector pMpGE_En01 (Tokushima University) |
| MN044709.1 | attR1, attR2 | Sterner, Z.R. et al., Genesis 57:e23329 (2019): Gateway destination vector pDXTP (University of Cincinnati) |
| LC122270.1 | attP1, attP2 | Kusano, H. & Shimada, H. (direct submission, Tokyo University of Science): Gateway donor vector pSSARL-GwBP |

## Gateway att sites

`attB1/2`, `attL1/2`, `attR1/2` and `attP1/2` are cut from the four academic vector records above,
at the places those records' own feature tables annotate, each on the plus strand so every entry
reads in the canonical orientation. The sites come to the expected lengths (attB 25 bp, attL 100 bp,
attR 125 bp, attP 232 bp) and each carries the core its number implies: `TTTGTACAAAAAAGC` for the
1 sites and `CTTTCTTGTACAAAG` for the 2 sites, the same cores `snappiesnap_bio::cloning` uses to
simulate BP and LR reactions (Hartley, J.L., Temple, G.F. & Brasch, M.A., Genome Res. 10:1788
(2000)).

The original Gateway vector records are company submissions, so they are not used (see "Which
records may be used" above). The sites are those of the published academic vectors listed, which
is also why the attP record is a laboratory donor vector rather than a catalogue pDONR.

## Fluorescent proteins (FPbase)

`fpbase_proteins.json` lists 976 fluorescent proteins from FPbase (https://www.fpbase.org),
built by `cargo xtask fpbase` on 2026-10-24 from the public API
(https://www.fpbase.org/api/proteins/). It keeps every protein whose sequence has only the 20
standard amino acids. Each entry holds the name, the protein sequence, the DOI of its primary
reference, its GenBank protein id (where FPbase gives one) and its FPbase page.

FPbase licenses its data under **CC BY-SA 4.0** (https://www.fpbase.org/terms/) and asks users
to credit the original authors of each protein. SnappieSnap does both:
- the file keeps that licence;
- a feature added from it carries the FPbase page, the licence and the primary reference
  (`doi:`) in its note.

FPbase also notes that some proteins may be covered by patents or other rights. Detecting
and naming a protein in a sequence grants no rights to use it.

## Boundaries from publications

Some elements are annotated in their records only as a point, or not at all. These take their
boundaries from a published paper (`citation` in the manifest), measured on the named record.
The generator checks that the cut reads exactly the sequence the paper gives (`check`), and a
feature added from such an element names the paper in a note.

## Downloaded vector library (not shipped)

`vector_records.json` lists only the accession numbers of 25 common vector records in GenBank.
Many of them were submitted by vector suppliers. On request, and after the user accepts its
terms (`feature_library_download` with `accept: true`), SnappieSnap downloads these records
from NCBI and builds a feature library from their annotations in the user's data folder.
The library is for personal use on that computer and may not be shared. Nothing from these
records is ever shipped with SnappieSnap.

## Adding a feature

1. Find a genome, academic or RefSeq record that annotates the element itself.
2. List the record under `records` in `manifest.json`, with why it may be used (who
   submitted it). Then add the feature with the record's key and location exactly, plus a word
   from its qualifiers as `expect`.
3. Or, when the record doesn't annotate the element: leave out `key`, give the `location` on
   the record, the `citation` (with PMID) and the published sequence as `check`.
4. Run `cargo xtask features`, then add a row above.
