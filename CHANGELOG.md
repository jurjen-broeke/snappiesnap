# Changelog

What each version of SnappieSnap added and fixed, newest first.

Versions are numbered 0.x.y. A version that adds features raises x. A version that only fixes things raises y. Each entry lists **Added** (new features) and **Fixed** (things that were broken or didn't work) separately.

## 0.41.1 (2026-11-26)

### Fixed

- **Feature names fit inside their arrows in the sequence view.** The feature bars were so thin that their names were squeezed against the edges. The bars are now taller, so a name has room inside its arrow.

## 0.41.0 (2026-11-26)

### Added

- **A navigator under the one-line sequence view.** When the sequence runs along one line, a strip along the bottom of the view shows the whole molecule: every feature as a bar in the same colour as in the sequence, your selection, and a box marking the stretch on screen. Click or drag in it to move the view there, and hover over it to see the position and the names of the features under the pointer. It shows the construct only; sequences aligned in the Align with pane are left out.

## 0.40.2 (2026-11-26)

### Fixed

- **The Align with pane no longer creeps up over the sequence.** Every time the mouse moved, the pane under the sequence grew a little taller, until it pushed the sequence view above it out of sight. It now keeps the height you give it.

## 0.40.1 (2026-11-26)

### Fixed

- **On Windows, `SNAPPIESNAP_DATA` and `SNAPPIESNAP_SETTINGS` were ignored when set to a path starting with `/`** (such as `/data/snappiesnap`), and SnappieSnap quietly used its usual folder in `%APPDATA%` instead. Such a path is rooted, not relative, and is now honoured as the documentation says.

## 0.40.0 (2026-11-25)

### Added

- **Check for updates now works.** The button in Preferences asks the project's release page for the newest version and says whether you are up to date, instead of reporting that no release has been published. Nothing is checked automatically, and nothing about you is ever sent.

## 0.39.0 (2026-11-25)

### Added

- **Align reads under the sequence.** The sequence view has an **Align with** panel (the new button under the layer toggles, or *Align reads under the sequence…* in the Align tab). List the reads and sequences to align there: add them with **Add files…**, pick one already open with **From a tab**, or drop files on the panel (traces dropped anywhere on the sequence view while it is open land there too). Reads stored in the `.dna` file are listed by themselves. Tick, untick and reorder them with **Up** and **Down**.
- **The whole alignment under the sequence**, in a pane you can make taller or shorter: the construct, then every ticked read in the list's order with its mismatches boxed, poor calls faint, its strand, its length, the part that aligned and its counts of mismatches and gaps.
- **Each read's trace under its own bases.** The chromatogram is stretched column by column, so every peak sits under the base it was called as — reverse reads drawn as the top strand reads them. Click a read's name to fold its trace away and keep one line of bases.
- **A navigator for the whole construct** at the bottom: where each read lies, its differences, and the stretch in view. Click or drag it to move there; swipe sideways or Shift-scroll in the pane to scroll along.
- Clicking a column selects that base in the sequence above, and a selection made above brings the alignment to it.

## 0.38.0 (2026-11-17)

### Fixed

- **A primer with changes designed into it is no longer cut short.** A repair or mutagenesis oligo carries the very differences it is meant to make, and the binding rule used to stop at the first one and report everything beyond it as a 5' tail that hangs free. A 100-base repair oligo carrying three changes showed only its last 49 bases. It now pairs along its whole length, across any number of single changed bases, as long as its 3' end — the part that decides priming — still pairs cleanly.
- A primer that pairs along its whole length is drawn as its real binding site rather than greyed out as a weak one. PCR is unchanged: polymerase still only extends from a site that pairs exactly, so no new products are invented.

### Added

- **A base that doesn't pair steps up out of the row**, with the row stepping around it, so a primer still reads as one piece and its changes stand out where they are.

## 0.37.0 (2026-11-17)

### Added

- **The sequence on one line.** The sequence view can run the whole molecule along a single line that scrolls sideways instead of wrapping it into rows (the ↔ button beside the other layer toggles). A ruler of positions runs along the top, and the row labels stay against the left edge as the residues scroll past. For following a feature, a primer or a read along the sequence, nothing wraps out from under you.
- **The read check's alignment has the same choice** — tick **One line** above it, and the construct and every read run along together.
- Only the residues in view are drawn, so one line of a 40 kb molecule costs no more than one line of a plasmid.
- **A selection that lands off screen brings itself into view.** Go to, Find, and commands an agent runs now scroll the sequence view to what they selected; a selection made with the pointer is on screen already, so the view stays where you put it.

## 0.36.0 (2026-11-16)

### Added

- **Reads can be tabs.** Check Reads lists the constructs and traces already open in tabs, to tick as reads beside any files. A trace dropped on the page or opened to look at its chromatogram can be checked straight away, and its own trace still shows under the differences it covers.
- **The Align tab puts the alignment beside the inputs** on a wide page, instead of squeezing it under them with the rest of the page empty — the same layout Check Reads got, for both To this sequence and Multiple.
- **Alignments use the width they are given** everywhere they are written out: a block is as many residues as there is room for rather than always 60, so there is less to scroll through on a wide window.

### Fixed

- **Add read files… did nothing in the web version**, which left no way to add a read at all and no way past the first step. It now opens the browser's file picker and takes several traces at once.

## 0.35.0 (2026-11-16)

### Added

- **Check Reads puts the alignment in a column of its own.** On a wide page the summary — how the reads cover the construct, how each one lies, and every difference — stays on the left, and the reads themselves are written out on the right, with the chromatograms of the chosen difference under them. The page's width is used instead of squeezing the alignment under everything else. A narrow window still stacks them.
- **The alignment uses the width it is given**, fitting as many columns to a block as there is room for (90 instead of 60 on a wide page), so fewer blocks to scroll through.

## 0.34.0 (2026-11-16)

### Added

- **Align a sequence that is already open.** Both alignments now take their sequences from the other open tabs: **To this sequence** has a **From a tab** choice beside From a file and Paste, and **Multiple** lists the other tabs to tick. No saving, copying or pasting to align two constructs you already have in front of you.
- **Load alignment files with a picker.** **Multiple** has an **Add files…** button that takes several files at once, and the pairwise alignment has **Choose…** for one. The files you load are listed beside any from the open folder.
- This matters most in the browser, where a page is given no folder to browse: until now the only way to align a file there was to open it, copy it and paste it in, one at a time.
- Agents get the same: `align_to_reference` takes a `document` (an open tab's id) in place of `file` or `residues`, and `align_sequences` takes `documents`.

## 0.33.0 (2026-11-16)

### Added

- **Check Reads lays the reads out under the construct**, column by column in blocks of 60, each read's bases written as they pair with the construct's top strand, so several reads of the same stretch can be read down the page. A base that differs from the construct stands out, and so does a gap where a read is missing one; extra bases in a read are marked on the join rather than pushing the other reads out of step.
- **A poor call is drawn faint.** Anything called below Q20 is dimmed, so the noise that survives trimming near a read's ends looks like noise and a well-called difference in the middle stands clear of it.
- **Everything points at the same difference.** Clicking a column of the alignment, a tick on the coverage map, or a row of the difference list picks that difference everywhere: the alignment moves to its block, the column is marked, and the chromatograms of every read covering it appear below.
- Only the stretch the reads actually lie on is drawn, so reads of one small region of a long construct don't leave screens of empty columns to scroll through.
- Agents get the same layout: `check_reads` takes `rows`, and each placed read then carries a `row` against the construct, a `row_quality` beside it, and `row_insertions` for the columns it has extra bases before.

## 0.32.0 (2026-11-16)

### Added

- **The product step shows the plan as well**, under the product's map, so the same picture of the reaction runs through every step: the vector it started from, the piece that went in, and the product as a circle whose black and red arcs are the backbone and the insert in their real proportion.
- **The plan warns you when something is open to doubt.** If both ends of the piece take the same cut it can go in either way round, so the plan says to expect both orientations and screen for the one you want. If the ends don't match the backbone at all, it says that before you run anything. And if the enzymes cut the vector into more than two pieces, it says the reaction keeps the largest, so you can check that is the backbone you meant.

## 0.31.0 (2026-11-16)

### Added

- **The plan under each cloning step now draws the pieces the step actually makes.** A plasmid is a circle; cutting it with two enzymes leaves a circle with a gap, with the small piece taken out drawn in that gap, between the two cut ends and to the same scale. A single cut simply breaks the circle open. A linear insert cut at both ends keeps its middle, with the flanks it loses shown greyed and set apart. The sizes come from the digest itself, so they are the pieces the reaction will really produce, and each one says how many base pairs it is.

## 0.30.0 (2026-11-16)

### Added

- **Point at an enzyme's name in the sequence view to see where it cuts.** Its cut through both strands is drawn as a thick line in a highlight colour, and its name turns the same colour, so it's easy to tell which cut belongs to which name when several sit close together. This works in document tabs and in the Cloning experiment's sequence view.

## 0.29.0 (2026-11-16)

### Added

- **Double-click an enzyme's name in the sequence view to add it to the cut list** in the Cloning experiment, just as on the map. In a document tab a double-click on an enzyme does nothing new.

### Fixed

- **A CDS's translation could be hard or impossible to read in the sequence view**, because the amino acids were drawn inside the feature's arrow, where its colour and its name covered them. Each translation now has a row of its own directly above its feature, in the text colour, and the features come straight under the two strands with the reading frames below them.

## 0.28.0 (2026-11-13)

### Added

- **Each cloning step shows its construct the way a document tab does.** The layer toggles are down the left edge, the Map / Sequence switch, the map's zoom (−, %, +, Fit, with pinch and drag) and the Auto / Circular / Linear choice are above it, and the sequence view shows the reading frames and genetic code you have chosen. The only difference is that it is read-only: you can look, select and copy, but a construct is changed in its own tab.
- **Pick from the construct instead of typing.** Double-click an enzyme to add it to that side's cut list, or a primer to put it in the forward or reverse field — by the strand it binds, since a pair is one primer on each strand, and choosing one sets the side to *as a PCR product*. Clicking a feature selects it, as in a document tab. The **Add…** list and the fields are still there for typing.

## 0.27.0 (2026-11-13)

### Added

- **Every enzyme that cuts the construct is on the map in the Cloning experiment**, not a shortlist, and the ones you have named for the step are drawn larger and in the accent colour. They also keep their labels when the map runs out of room, so the enzymes the step is about can't be crowded out by the others.
- **Double-click an enzyme on the map to add it to the cut list.** The **Add…** list is still there; this is the quicker way to pick from what the construct actually offers. Naming the same enzyme twice doesn't double it up.

### Fixed

- **The Back and Continue buttons stayed below the window** on a large display: the two columns took the full height and pushed the step's buttons off the bottom. They now keep their room.
- **ORFs and primer binding sites didn't show in the Cloning experiment.** The toggles were there from 0.26.0 but nothing was behind them: the pane passed empty lists. It now finds both for the construct it is showing.

## 0.26.1 (2026-11-13)

### Fixed

- **Selecting and copying work in the Cloning experiment's sequence view.** The pane is read-only, but it was refusing everything: a selection is now kept (per construct, since a range means nothing in another molecule) and Copy, Copy as FASTA and Copy reverse complement put the same text on the clipboard that the Copy command does. Only changes — typing, cut, paste — are refused, and the note now says what does work.
- **The plan shows only the step you are on.** It was drawing the whole reaction at every step; that picture belongs to the last one. Step 2 shows the vector and what cutting it gives, step 3 the piece that goes in and how it is prepared, and the product step the whole reaction as before. Row names widen their column instead of running over the drawing.

## 0.26.0 (2026-11-13)

### Added

- **The Cloning experiment shows a construct the way a document tab does.** Every step now has **Map** / **Sequence** buttons and the same layer toggles — features, restriction sites, primer binding sites, ORFs and the ruler — so you can read the sequence of the vector or the insert, with its features and cut sites, without leaving the experiment. The toggles are the ones the document views use, so a layer means the same thing wherever you switch it on.
- **The product step gets the same treatment**: the full width of the window, the view switch, the layers and the enzymes that cut the product, instead of the small fixed map it had.
- **The plan says whether each piece is circular or linear.** The vector, what goes into it and the product each carry the word at the end of their row, and a circular piece is drawn with its ends joined. A construct that hasn't been read yet says nothing either way rather than guessing.

## 0.25.0 (2026-11-13)

### Added

- **The Cloning experiment is laid out like the gel experiment.** The choices sit in a column on the left (a third of the window) and the construct fills the rest, so the map is big enough to work from. **The enzymes are drawn on it**: the ones you typed in the fields, or — before you have typed any — the unique cutters you are most likely to own, up to twelve of them, so the map helps you choose instead of burying you in labels.
- **A plan under the map** shows what the reaction will do: the vector as a black line, everything going into it in red, drawn to scale from the real lengths, with the ends named by the enzymes that make them and the product closing into a circle. It follows the technique: Gateway names att sites rather than enzymes, mutagenesis has no insert row, and an assembly's extra parts each get their own piece. A construct that hasn't been read yet is drawn as an outline rather than at a made-up size.

## 0.24.0 (2026-11-13)

### Added

- **Build your own enzyme set inside the app.** Tick enzymes in the Enzymes tab, give the set a name and press **Save as set** — no file, no dialog, so it works in the **browser** as well as on the desktop. **Tick all** fills the ticks from whatever the table is showing, so an existing set is edited by listing it, changing the ticks and saving it under the same name (the button says **Replace set**, and names the set it would replace). The new set becomes the one the map and the sequence view draw from.
- **The browser can import an enzyme set file too**, through its own file picker (**Import a set…** in the Enzymes tab). Sets made or imported in the browser last for the visit — a web page cannot write files that outlive it — and the tab says so.
- An enzyme your set names that the database doesn't have is kept in the set and reported rather than dropped, so a lab that owns a discontinued or unusual enzyme can still name it; downloading REBASE makes it usable.

### Fixed

- The Enzymes tab no longer says "587 enzymes (587 commercially available)": the built-in set is every enzyme a supplier sells, so it now reads "all sold by a supplier". The count still appears for REBASE, which holds research enzymes nobody sells.

## 0.23.0 (2026-11-13)

### Added

- **A cloning experiment you can run with nothing of your own open.** The Cloning experiment's first step has a **Try it with the sample constructs** button: it opens pSample and demoInsert and fills in a restriction clone with EcoRI and BamHI, so you can step through the vector, the insert and the product and open a real construct. This is the whole of a first try in the **web version**, where a browser page has no files to choose from.
- **pSample has a polylinker.** The sample plasmid gained a 38 bp polylinker with six enzymes that cut it exactly once — EcoRI, BamHI, SalI, XhoI, NotI and PstI — written over a stretch that held nothing, so the plasmid is still 3,200 bp and every feature is where it was. The new **demoInsert** is a 762 bp fragment with EcoRI and BamHI ends that drops into it, giving a 3,950 bp product with the insert in the polylinker and every original feature carried over.
- `sample_open` opens the samples for agents and scripts too (`which`: `plasmid`, `insert` or `cloning`). Like every other command, it refuses to replace unsaved work unless asked.

## 0.22.0 (2026-11-13)

### Added

- **The built-in enzymes now cover everything a supplier sells.** The set grew from 260 to 587: the rule was "sold by at least two suppliers", which quietly dropped the single-supplier catalogue names labs actually buy, so an imported enzyme set lost them unless REBASE had been downloaded. Of the 84 enzymes in a real lab's set, 83 now resolve straight away (the one left, BlnI, is a name REBASE lists with no supplier at all; AvrII and XmaJI are the sold names for its site). Names such as MscI, HindII, MunI, KspI, AflIII, BsiWI, SexAI, XcmI and Hpy188III are built in.
- **One name per site on maps.** Where several enzymes cut a site the same way, the map and the sequence view now label it with a single name — the one your own enzyme set uses, otherwise the one most suppliers sell — and the hover card lists the other names. Labels stay short as the enzyme set grows, and they speak your lab's vocabulary.

## 0.21.0 (2026-11-13)

### Added

- **Clone from the constructs you have open.** Every step of the Cloning experiment can take a construct straight from an open tab, listed under **Open now** in its picker above the files of the open folder. This is what the **web version** needs: a browser gives the page no file paths, so until now the Cloning experiment had nothing to choose from there and could not be used at all. It also means a construct you have just edited, or a product of an earlier reaction, is cloned as it now stands, without saving it first; edit it in its tab and the experiment follows. The commands do the same for agents: `restriction_clone`, `gateway_clone`, `assemble`, `topo_clone` and `mutagenesis` take `vector_document`, `insert_document`, `acceptor_document`, `donor_document`, `template_document` or a part's `document` — an open document's id, from `documents` — instead of a file.

- **The cloning experiment offers the enzymes you have.** Every enzyme field in the Cloning experiment (the vector's and insert's enzymes for restriction cloning, the enzymes a part is cut with for Gibson and In-Fusion, and Golden Gate's Type IIS enzyme) now has an **Add…** list of the enzymes that cut *that* construct, fewest cuts first, with the number of sites. When an enzyme set of your own is chosen, the list holds only what the set has, so the reaction you design is one you can run; without a set, every commercial enzyme is offered and the field says so. Typing is never restricted: an enzyme the set hasn't got is accepted and marked *you'd have to buy or borrow it*, and a name the enzyme database doesn't know is reported before the reaction is run, with a pointer to downloading REBASE.

## 0.20.0 (2026-11-06)

### Added

- **Diagnostic Digest searches the enzymes you have.** An enzyme set of your own is an inventory of the fridge, so when one is chosen the designer searches only it, and every digest it offers is one you can run today. Tick **Include enzymes I haven't got** to search everything: those suggestions are marked *not in your set* and listed after the ones you can run. Enzymes in your set that the enzyme database doesn't know are named, with a pointer to downloading REBASE, instead of quietly dropping out of the search. `diagnostic_digest` takes `set` and `include_out_of_stock`.

## 0.19.0 (2026-11-06)

### Added

- The oligo library imports a whole folder at once. A file holding nothing but a sequence is named after the file, which is how a lab keeps its cloning primers: one file per primer, the name in the file name. Lists and FASTA files in the folder are read too, and anything else there is left alone. The Lab Oligos section has an **Import folder…** button beside **Import oligo list…**, and the status line says how many files were read and how many held no oligos.

## 0.18.1 (2026-11-05)

### Fixed

- Gateway no longer simulates reactions that cannot happen. attB, attL, attR and attP all share the same core, so until now any two constructs carrying att sites could be "recombined" — an entry clone with an expression clone, say — and the result looked convincing. Which kind each site is, is now read from the arms around its core: LR needs an entry clone (attL) and a destination vector (attR), BP an attB insert and a donor vector (attP), and anything else is refused with a message naming what each construct carries and what the reaction needs. Handing the two over the wrong way round says so. Sites whose arms are too unusual to name are still simulated, with a warning that the reaction wasn't checked.

## 0.18.0 (2026-11-05)

### Added

- **Cloning is now an experiment of its own**, in the Experiments menu, built in four steps: the technique (restriction, Gateway LR and BP, Gibson/HiFi, In-Fusion, Golden Gate, TA and TOPO, site-directed mutagenesis), the vector, the insert, and the product. Both constructs are picked from the folder you have open (and its subfolders) or from disk, and each is shown as a map before you go on, so no construct has to be open first. Step 3 runs the reaction as you set it up and says whether the two fit, or what is wrong, before you commit to it. Step 4 shows the product and opens it as a new construct with its construction history. An assembly can take more than two parts, and mutagenesis takes its two primers in step 3 instead of a second construct. This replaces the Cloning tab, which is gone.
- Site-directed mutagenesis can work on a plasmid in a file (`template`), not only on the open document.
- A new `file_info` command reports what a sequence file holds — name, length, topology, how many features and primers, and the molecule itself for a map — without opening it as a document.

## 0.17.0 (2026-11-05)

### Added

- Simulate Gel is now one item in the Experiments menu, covering both kinds of gel. It opens on a new first step that asks what you want to simulate: a gel of several constructs, or a diagnostic digest. The steps after it follow your choice, and you can go back to it from the constructs step without losing the constructs you picked. Diagnostic Digest no longer needs a menu item of its own; Check Reads keeps one, since it runs no gel.

## 0.16.0 (2026-11-05)

### Added

- The eight Gateway recombination sites are detected out of the box: attB1 and attB2, attL1 and attL2, attR1 and attR2, attP1 and attP2. Detect Features now names them in entry clones, donor and destination vectors, and expression clones, with no feature library of your own. Each site is cut from an academic laboratory's vector record, at the place that record's own annotation gives, and an added feature names the record it came from.

## 0.15.1 (2026-11-05)

### Fixed

- In the browser, Preferences → Files and tabs no longer offers a folder panel switch that did nothing. It explains that a web page can't see your folders.

## 0.15.0 (2026-11-05)

### Added

- SnappieSnap runs in the browser. The web version is the same app compiled to WebAssembly: a folder of static files that any web server can host, with nothing to install for users (see docs/web.md).
- Open files with the browser's file picker, or drop them on the page. Everything runs on your computer; files are never uploaded.
- Save downloads the file. Save As, Save Copy without History, experiments and map export ask for a file name and then download it.
- The browser asks before you close or reload the page while something is unsaved.

## 0.14.0 (2026-11-02)

### Added

- The Sequence view writes each primer out base by base against the strand it pairs with: forward primers 5'→3' under the top strand, reverse primers 3'→5' under the bottom one, so every base sits over its partner. Mismatched bases are shown in red, bases that don't pair in lowercase grey, and overlapping primers stack in rows of their own.
- A 5' tail (a restriction site or overlap added to a primer) is shown as a small loop with its length, or written out in full: Preferences > Primers > Sequence view.
- For agents: `ui.set` takes `primer_tails` (`simplified` or `full`).

## 0.13.0 (2026-11-02)

### Added

- Feature detection now also finds copies with inserted or missing bases, not only substituted ones, within the same allowance (2% of a feature's length by default). The Features tab says how each near match differs, and a feature added from it gets a note such as "3 of 861 bases differing (2 inserted, 1 missing)".
- For agents: `detect_features` and `add_detected_features` take `indels` (default true; false for substitutions only) and list `substitutions`, `insertions` and `deletions` for each find.

## 0.12.0 (2026-11-02)

### Added

- Export maps: the Export… button in the Map view saves the whole map as SVG, vector graphics for papers and posters, or as a PNG image at twice the size. The map is saved exactly as shown (layers, enzymes, colours and background) but unzoomed and without the selection; in SVG files the labels stay editable text. Circular maps are saved square, linear maps wide.
- `snappiesnap-cli render file.dna map.svg` saves the map as SVG without opening a window.
- For agents: `ui.export_map`.

## 0.11.0 (2026-11-02)

### Added

- Sequencing reads matched to a construct in SnapGene are now read from its .dna file: traces and plain sequences such as whole-plasmid assemblies. Opening such a file says how many reads it has.
- Check Reads can include these reads ("Include the reads stored with …"); they are used as trimmed in SnapGene. Every read in the results has an Open button that opens it in a tab of its own, traces with their chromatogram in the Trace tab.
- The differences of a read check can be limited to those seen in calls of quality 20 or better, which hides most of the noise at the ends of reads.
- For agents: `matched_reads`, `open_matched_read`, and `check_reads` with `matched`.

## 0.10.1 (2026-11-02)

### Fixed

- The construction history of some .dna files couldn't be shown ("the two lengths in the sequence header differ"). These files store earlier steps that contain unknown bases (N) in a form SnappieSnap didn't read yet; such histories now open, with the N bases in place. Histories that SnappieSnap writes can now hold N bases too.
- Steps that SnapGene keeps hidden in a history, without a sequence of their own, no longer stop the rest of the history from being read.

## 0.10.0 (2026-11-02)

### Added

- Primer sites that match only at the 5' end: where a primer's first bases (15 by default) pair with the sequence but its 3' end doesn't, the primer can anneal but can't be extended. These sites now show on the map and in the sequence view, and their hover card draws the paired 5' part and the unpaired 3' end. Turn them off, or change the number of bases, in Preferences › Primers.
- For agents: `primer_binding` takes `five_prime_on` and `five_prime_matched`, and `primer_sites` marks such sites `five_prime_only`.

## 0.9.0 (2026-11-02)

### Added

- Check Reads (Experiments menu): choose a construct and its sequencing reads (.ab1, .scf, FASTA or GenBank files, or every trace in the folder at once). Each read is trimmed of its poor ends and placed on the construct, on either strand. The results show a map of what the reads cover, with the stretches no read covers in red, and a table of the reads with their strand, place, identity and the part kept.
- Every difference from the construct is listed once for all reads: confirmed (every read covering it shows it, and at least two do), disputed (other reads show the construct's base) or seen in one read, with the quality of each read's call. Click a difference, or its mark on the map, to see the trace of every read covering it, drawn as the construct's top strand reads. The first disputed difference opens by itself.
- Read checks are saved and reopened like other experiments (.sslab); the read files are kept as links, not copies.
- For agents: `check_reads`.

## 0.8.0 (2026-11-02)

### Added

- Diagnostic Digest (Experiments menu): choose the constructs to tell apart, such as a clone and its empty vector, and SnappieSnap suggests restriction enzymes, alone or in pairs, whose digests give clearly different bands on a gel. With one construct it suggests digests whose pattern confirms it. Each suggestion lists the fragment sizes of every construct and what to watch out for (faint small fragments, a buffer for double digests). Use digests every lane and shows the simulated gel.
- Enzymes that cut the same way are listed once, under the most specific recognition site and the name your lab uses (from your enzyme sets) or most suppliers sell; the others are shown on hover.
- Saved experiments remember that they are a diagnostic digest.
- For agents: `diagnostic_digest`.

## 0.7.0 (2026-11-02)

### Added

- Go to a position or range: press Cmd+L (Ctrl+L on Windows and Linux) in the Map or Sequence view, or click Go to… in the status bar, and type a number such as 1200 or a range such as 100-250. On a circular molecule a range like 4000-200 runs across the origin.
- The status bar describes the selection, or the whole molecule when nothing is selected: GC content, the Tm of stretches of 8 to 60 bases (worked out as for primers), and the molecular weight of one strand or of both. Proteins show their mass.
- For agents: `go_to` and `sequence_info`.

### Fixed

- The status bar of a protein gave its length and selection in bp; it now says aa.

## 0.6.1 (2026-11-02)

### Fixed

- Experiment files now end in .sslab ("SnappieSnap lab") instead of .ssexp. The format is unchanged, and .ssexp files saved by 0.6.0 still open.

## 0.6.0 (2026-11-02)

### Added

- Save experiments: an experiment tab has Save and Save As… buttons that write a SnappieSnap experiment file (.ssexp). The file holds the experiment's name, gel settings and lanes with their methods. Each construct is stored as a link to its file and as a copy, so the experiment can be moved or shared.
- Open experiments with Experiments › Open Experiment…, the Open button, the folder panel or by dropping a .ssexp file on the window. Each construct opens from its own file when that file is found and unchanged; otherwise the copy saved in the experiment opens, and the status bar says why. Constructs that are already open are reused.
- The status bar of an experiment tab shows whether it is saved or has unsaved changes.
- The .ssexp format is documented in docs/experiment-format.md, so other software can read and write it.
- For agents: `experiment_save` and `experiment_open`, and in the app `ui.open_experiment` and `ui.save_experiment`.

### Fixed

- The status bar of an experiment tab counted lanes as constructs; it now counts each construct once.

## 0.5.0 (2026-11-02)

### Added

- Sequencing reads: open .ab1 (Applied Biosystems) and .scf files. The read opens in a new Trace tab that shows the four dye signals, each base above its peak, and a bar for each base's quality. Low-quality bases are dimmed, and the poor ends of the read are shaded. Width and Height buttons zoom the trace.
- Select bases by clicking or dragging in the Trace tab; the selection is shared with the Sequence view. Bases you type in are marked with a dot because no peak backs them. Reverse-complementing a read flips its trace.
- Trim ends (Trace tab, or Actions › Trim Low-Quality Ends) deletes the shaded low-quality ends of a read in one undoable step.
- Align to This Sequence trims a read's low-quality ends before aligning it, so a read file can be checked against its construct directly.
- Saving an edited read asks for a new file (.dna, GenBank or FASTA). These keep the sequence but not the chromatogram, and the app says so. An unchanged read saves back to its original file exactly.
- For agents: `trace_info` (read quality, the stretch trimming keeps, and each base's call and signal) and `trace_trim`.

## 0.4.3 (2026-10-30)

### Fixed

- Circular maps: labels no longer run off the top of the map when one side has many of them; in the sample plasmid, "marker" and the enzymes near it were drawn outside the window. Labels now always stay within the map. When there are more than fit, some enzyme labels are left out (feature labels last), and a line at the bottom says how many.
- Circular maps: long labels are no longer cut off at the left and right edges. The map leaves room beside the circle, and a label that still doesn't fit is shortened with "…".

## 0.4.2 (2026-10-30)

### Fixed

- Primer binding rules now depend on the length of the sequence the primer binds to, not on the primer's length: one set of rules for short sequences and one for long sequences such as chromosomes (long from 1,000,000 bp by default). Defaults: short sequences 10 matched bases with an isolated mismatch allowed; long sequences 15 matched bases with no mismatch; at least 40 °C for both.
- "Allow isolated mismatches" now allows one mismatch, not one per 10 bases.
- Rules saved by 0.4.0 or 0.4.1 go back to the defaults once.

## 0.4.1 (2026-10-30)

### Fixed

- Primer binding rules: short and long primers now each have their own matched bases, isolated-mismatch setting and minimum Tm, instead of sharing the mismatch and Tm settings. Rules saved by 0.4.0 go back to the defaults once.

## 0.4.0 (2026-10-30)

### Added

**Primer binding rules (Preferences > Primers)**
- Set how many bases at a primer's 3' end must pair with the template, separately for short and long primers, and from what length a primer counts as long (defaults: 10 bases, and 15 bases from 25 bases on).
- Allow single isolated mismatches: one mismatched base with pairing bases on both sides, at most one per 10 bases and never in the last 3 bases (on by default).
- A minimum Tm for binding sites, on at 40 °C by default. Weak sites whose paired stretch melts lower are no longer shown.
- The rules apply to the map and sequence view, the Primers tab, the lab oligo search, and PCR when it also uses weaker sites. Exact sites over 15 bases or more (or the whole primer) always count. The rules are remembered between sessions.

## 0.3.0 (2026-10-30)

### Added

**Preferences**
- Files and tabs: show the folder panel (moved here from General) and choose its folder.
- DNA: what the map and sequence view show (features, restriction sites, primer binding sites, cut positions, ruler, reading frames), what the ORF layer shows, and the genetic code. The reading-frame choice and the genetic code are now remembered between sessions.
- Enzymes: which enzymes the map and sequence view show (all or one of your enzyme sets, how often they cut, site length, commercially available only), which enzymes the Enzymes tab lists, the enzyme database (download REBASE or go back to the built-in enzymes), and your enzyme sets (import or remove).
- General, Updates: a "Check for updates" button. SnappieSnap only looks for a new version when you press it. Until releases are published, it says there is nothing to check against yet.

## 0.2.0 (2026-10-30)

### Added

**Preferences**
- A new Preferences window: categories in a list on the left (General, Files and tabs, Export, DNA, RNA, Gels, Enzymes, Primers) and the chosen category's settings on the right, in groups. Preferences reopens on the category you last used.
- General, Window: the interface size (75 % to 200 %, also with Ctrl/Cmd + and -), showing the folder panel, and reopening the window at the size it had last time (on by default).
- General, Colour theme and background: the theme (like the system, light or dark) and the Map and Sequence background for each theme.
- General, Updates: shows the version you are running. Checking for new versions comes later.
- The other categories are in place, and their settings will follow in later versions.

## 0.1.0 (2026-10-29)

The first numbered version. It covers everything built so far.

### Added

**Version number**
- The window title shows the version, for example "pSample — SnappieSnap 0.1.0".

**Files and documents**
- Open and save SnapGene `.dna` files. A file opened and saved without changes stays byte for byte identical, and edits change only what they touch.
- Open, save and convert GenBank, EMBL, GFF3, SBOL 2 and FASTA files. GFF3 annotations can be added onto an open sequence.
- Several documents open at once, in tabs; each tab remembers how it was being viewed.
- A folder panel to browse, open, rename, move and copy sequence files.
- Save with or without the construction history.
- All of SnappieSnap's own files (settings, libraries, downloads) are kept in a folder of the user's own.

**Editing**
- Type, select, copy, cut and paste in the sequence view, also across the origin of circular molecules.
- Find sequences (both strands, IUPAC codes, up to 3 mismatches) and feature or primer names.
- Undo and redo with a named history; reverse complement.

**Maps and the sequence view**
- Circular and linear maps with features, restriction sites, primers, ORFs and a ruler, each switchable.
- Map zoom, and two documents side by side.
- The sequence view shows both strands, features, restriction sites, CDS translations and reading-frame translations.
- Preferences (gear in the toolbar): a background colour for the Map and Sequence views, one for each theme.
- A colour picker to give a feature its own colour. A feature's stored colour, set in SnappieSnap or read from a `.dna` file, is what the map and views show.

**Features**
- A Features tab with a sortable, searchable table and a feature editor.
- Feature detection: from a library learned from your own files, a built-in set from public INSDC/RefSeq records, fluorescent proteins from FPbase, and near matches that differ by a few bases.
- A vector feature library that can be downloaded from NCBI on request.

**Restriction enzymes**
- A built-in set of 260 commercially available enzymes, and the full REBASE database downloaded on request.
- Choose which enzymes the map and sequence view show: unique, 1–2 or 1–3 cutters, minimum site length, sold only, or a list of names.
- Lab enzyme sets: import the lab's own enzyme list (as exported from SnapGene) and show only those enzymes on the map and in the sequence view. The Enzymes tab lists the lab set by default, and the choice is remembered.
- An Enzymes tab with every enzyme's cuts and digests of the ticked ones.

**Primers and PCR**
- Primers from `.dna` files and where they bind; primers and notes can be edited and saved into `.dna` files.
- Melting temperatures, primer design for a selected region, and hairpin and dimer checks.
- PCR simulation; a product opens as a new document.
- A lab oligo library: import the lab's oligo list once, see which oligos bind the open construct, and add the ones you tick as primers.

**Cloning and construction history**
- Restriction cloning, Gateway, Gibson/In-Fusion, Golden Gate, site-directed mutagenesis, and TA/TOPO cloning.
- The construction history of `.dna` files is shown as a list or a diagram. Products made in SnappieSnap record their own history, which is written into `.dna` files.

**Gels and experiments**
- Simulated agarose gels for restriction digests and PCR products, with ladders.
- Experiments (toolbar): Simulate Gel puts several constructs side by side on one gel, each lane with its own method (uncut, digest or PCR). Each experiment opens in a tab of its own.
- Save a gel as a PNG or JPEG image, as shown or inverted (dark bands on white), with a legend naming every lane.

**Alignment and BLAST**
- DNA and protein alignment (BLOSUM62), alignment of sequences to the open construct, and multiple sequence alignment.
- BLAST a selection at NCBI in the browser.

**Agent control**
- An MCP server and a JSON control channel: every command, the document, the window state, offscreen renders and gel images are available to agents and scripts. Sessions can be recorded and replayed.

**Desktop**
- An application icon for Windows, macOS, Linux and the web.
- Drafts of the installers: Linux `.deb` and `.rpm`, Windows `.exe` and `.msi`, and a macOS `.dmg`.

### Fixed

- A feature's colour no longer changes between a construct, its PCR products and its history steps. It used to depend on the feature's position in the list.
- Features made from `.dna` inputs keep their colour when saved. They used to lose it, and showed grey in SnapGene.
- Undo uses far less memory on large sequences: 100 edits on a 10 Mb sequence went from about 1 GB to 70 MB.
- Find: a mixed-case query such as "attB" now searches names only, instead of also matching hundreds of sequence hits.
- GenBank import gives non-coding features (promoters, origins and the like) the same direction SnapGene does, and reads the `/direction` qualifier.
- GenBank import hides only the record-level `source` feature, not other `source` annotations. Editing a file keeps the spacing of its LOCUS line.
- `.dna` files: line breaks and tabs inside feature values are saved correctly, and saving a document with a description (for example one converted from GenBank) no longer fails.
- EMBL: a record without a sequence block no longer gains one when saved.
