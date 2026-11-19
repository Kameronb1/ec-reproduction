# Enterochromaffin cell analysis: reproducing the figure panels

This repository recreates the figure panels of the paper, main and supplementary, from Figure 1C on, starting from deposited input data. For each panel it runs four steps:

1. The analysis script computes the panel's numbers.
2. Those tables are compared with the tables the published panel was drawn from.
3. The plot script draws the panel from those checked tables.
4. The drawn image is compared with the published panel.

It covers 43 panels and 5 keys that several panels share: Figures 1C to 1E, 2A to 2E, 3A to 3C, 4A to 4E and 5A to 5C, and Supplementary Figures S1 to S9. Figures 1A and 1B are not reproduced. The supplementary tables are not reproduced either.

## The studies

The integrative analysis uses 15 public single-cell studies, with the cell types their authors published (the published labels):

- primary intestine: Elmentaite 2020 + 2021 (one study), Yu 2021a, Wang 2020, Burclaff 2022, Hickey 2023 and Egozi 2021;
- human intestinal organoid (HIO) grafts: McCauley 2023;
- stem cell-derived islets (SC-islets): Augsornworawat 2020, Balboa 2022, Du 2022, Augsornworawat 2023, Zhu 2023 and Siehler 2024;
- adult primary islets: Fasolino 2022 and Bandesh 2026.

The seven pediatric donors of Elmentaite 2020 with Crohn's disease are left out of every analysis, including the fit of the four-marker EC classifier; `shared/gut_reference_roster/gut_reference_roster.json` in the deposit lists them.

A few panels also use other public data: Veres 2019 (SC-islets; Figures S7 and S8A), the adult human brain atlas of Siletti 2023 (serotonin neurons; Figure 2B), two fetal pancreas datasets, Olaniru 2023 and Yu 2021b (Figure S6 only), and published gene lists from two preprints, Liu 2025 and Fiancette 2026 (Figure 3C). `upstream/SOURCES.tsv` lists every source with its accession and citation.

The label-free four-marker EC classifier is fitted inside the workflow (Figures S2 and S3). Its fitted constants are a result, not an input, and Figures 1E and 4B use them.

## What you need

- The input data deposit: Zenodo DOI [10.5281/zenodo.24093050](https://doi.org/10.5281/zenodo.24093050), 245 files, 4.79 GB unpacked (the zip is 2.14 GB). Unpack it into a folder.
- Three inputs are not in the deposit and are downloaded from their original source instead: MSigDB 2026.1 (Zenodo record 18968178) and two supplementary gene lists of bioRxiv preprints (Figure 3C). Their URLs and checksums are in `inputs/` (rows with origin `url`).
- Python 3.13 and [Snakemake](https://snakemake.readthedocs.io) 9.23.1: `conda env create -f envs/python.yaml`.
- R 4.6.1 with the packages in `envs/r_packages.tsv` (CRAN plus Bioconductor 3.23): `Rscript envs/install_r_packages.R`, then `Rscript envs/check_versions.R`. R must have the cairo graphics device (`capabilities("cairo")` is `TRUE`), and the Arial font must be installed.
- To install the exact version of every package, dependencies included: `pip install -r envs/python_lock.txt` for Python, and `renv::restore(lockfile = "envs/renv.lock")` for R (needs the [renv](https://rstudio.github.io/renv/) package).
- Memory: 8 GB is enough. Several steps run at once, so on a computer with little free memory use fewer cores (for example `--cores 2`).

## Run it

From the repository root:

```bash
python workflow/scripts/common/fetch_url_inputs.py /path/to/unpacked/deposit
snakemake --cores 4 --resources heavy_r=1 --config data_dir=/path/to/unpacked/deposit
```

The first command downloads the three inputs that come from their original source into the same folder and checks each one's SHA-256. If a server refuses the download, the script says so: download the file by hand from the URL in `inputs/` and put it at the path given there. `heavy_r=1` makes the memory-heavy R steps run one at a time. A full run of all figures took about 52 minutes on one core (`--cores 1`).

On Windows, keep the repository and the results folder in a short path (for example `C:\ec\`). Snakemake keeps a record file for every output, named after the output's full path, and a long path can exceed the 260-character limit of Windows; the run then stops with "Error recording metadata for finished job".

Snakemake works through these steps:

1. Check every input file against the checksum in `inputs/<figure>.tsv`. It stops if anything is missing or changed.
2. Run the analysis scripts, which write tables to `results/<figure>/`.
3. Compare each table with the published one in `expected/<figure>/`. A difference stops the run, and the log shows the first differences. Each passing check is recorded in `results/<figure>/checks/`.
4. Draw each panel as a 600 dpi PNG on a white background: `results/<figure>/panels/<panel>.png`. A panel is drawn only after every table it reads has passed its check.
5. Compare each drawn panel with the published one in `expected/<figure>/panels/`. A difference stops the run, and the log shows how large it is. Each passing check is recorded in `results/<figure>/checks/panels/<panel>.png.json`.

To build one figure only, run `snakemake --cores 4 --resources heavy_r=1 --config 'figures=["figure3"]' data_dir=...`. Some figures use another figure's results, and building them alone runs that figure's analysis too:

- Figures 1E and 4B use the classifier fitted for Figures S2 and S3 (`figure_s2_s3`).
- Figures 3, 5, S4 and S6 use the gene groups of Figure 2D (`figure2`).
- Figure S1A draws a table made by the Figure 1 rules, Figures S7 and S8B-C draw tables made by the Figure 4 rules, and Figure S9 draws tables made by the Figure 5 rules. Figure S8A also reads the sample rows of Figure 4, and Figure S4 uses the color range of Figure 3B.

## Find the code for a panel

`PANELS.tsv` has one row per panel. Each row follows the panel from its inputs to its image:

- `shows`: what the panel shows.
- `inputs`: the input manifest.
- `rules`: the Snakemake rules file.
- `scripts`: the analysis scripts.
- `tables`: the tables in `results/` that the panel is drawn from.
- `checked_against_published`: how many of those tables have a published copy to check against ("no table" for a key that draws no data).
- `plot`: the plot script.
- `image`: the drawn panel.
- `published_image`: the published panel.

For example, to see how Figure 2C is made, find its row. The analysis script it names computes the numbers, and `workflow/scripts/figure2/plot/2C.R` draws them.

`PANELS.tsv` is built from `panels/<figure>.tsv`. After editing those, rebuild it with `python workflow/scripts/common/build_panels_index.py`.

Each `<figure>` folder holds one figure:

| `<figure>` | Panels |
|---|---|
| `figure1` | 1C, 1D, 1E |
| `figure2` | 2A, 2B, 2C, 2D, 2E, and the keys of 2C and 2D |
| `figure3` | 3A, 3B, 3C |
| `figure4` | 4A, 4B, 4C, 4D, 4E, and the key of 4D-E |
| `figure5` | 5A, 5B, 5C |
| `figure_s1` | S1A, S1B, S1C, S1D |
| `figure_s2_s3` | S2A, S2B, S2C, S2D, S3A, S3B, S3C |
| `figure_s4` | S4A, S4B |
| `figure_s5` | S5A, S5B, and their key |
| `figure_s6` | S6 |
| `figure_s7` | S7 |
| `figure_s8` | S8A, S8B, S8C, and the key of S8B-C |
| `figure_s9` | S9A, S9B, S9C, S9D |

Figures S7 and S9 have no analysis rules of their own; their tables come from Figures 4 and 5. Panel names, plot scripts and panel images carry the paper's figure numbers and panel letters: `S9A.R` draws Figure S9A. Short prefixes in other file names, such as `c5_` or `s2_`, only group files and are not figure numbers.

| Folder | Contents |
|---|---|
| `workflow/scripts/<figure>/` | the analysis scripts |
| `workflow/scripts/<figure>/plot/` | one plot script per panel, `helpers.R` with the figure's shared drawing code, and small drawing files (hand-placed label positions, row layouts, the terms chosen for display) |
| `workflow/scripts/common/` | the input check and download, the table and image comparisons, the panel-drawing driver and the `PANELS.tsv` builder |
| `workflow/rules/<figure>.smk` | which inputs each analysis script reads and which tables it writes |
| `inputs/<figure>.tsv` | every input file: checksum, size, what it is, how it was made |
| `panels/<figure>.tsv` | each panel's tables, scripts, plot script, image name and place on the figure page |
| `expected/<figure>/` | the published tables: the tables the published panels were drawn from |
| `expected/<figure>/panels/` | the published panels |
| `config/` | run settings, and the study name for each study code used in tables and file names |
| `envs/` | the Python environment, the R package versions, scripts to install and check them, and lock files with the exact version of every package |
| `upstream/` | how the deposited inputs were made from raw public data, and the list of sources |

The folders of the data deposit (`shared/`, `figure1/` to `figure5/`, `figure3_supp/`, `figure_s6/`) group the input files; a file's folder does not say which figure reads it. `inputs/<figure>.tsv` lists the inputs of each figure.

## How a panel is drawn

Each plot script defines `draw_panel()`. It reads only two things: tables in `results/` (those of its own figure and of the figures named above) and the drawing files in its own folder. The plotting code shows how each panel is drawn: filtering, aggregation, smoothing, error bars and axis transformations.

`workflow/scripts/common/draw_panel.R` opens a device the size of the figure page (7.5 × 10 in at 600 dpi, cairo). It draws the panel at its place on that page and cuts out the panel's box, given in `panels/<figure>.tsv` as x, y, width and height in pixels.

Drawing on the whole page keeps every text and point position the same as in the published figure. A few panels follow how their published pages were made instead:

- Figures 3A to 3C and S6 are drawn on a canvas of their own and then placed on the page.
- The panels of Figures 1 and S1 are first drawn once on a 72 dpi page, as their published pages were, because ggplot2 caches font measurements per device type and the text positions depend on that cache.
- Each panel of Figure S9 lays out the whole page, because each panel's place depends on the heights of the panels above it, and then draws only its own panel.

## How closely the images match

The `image_check` setting in `config/config.yaml` picks the comparison:

- **`exact`.** Every pixel must be identical. This holds on the setup the figures were made with: R 4.6.1 on Windows, the cairo PNG device, Arial, and the package versions in `envs/r_packages.tsv`. On that setup all 48 drawn images (43 panels and 5 keys) match exactly.
- **`tolerance`** (the default). Other systems draw text and line edges slightly differently, so pixels shift without the picture changing. Both images are averaged over 24 × 24 pixel blocks (1/25 in). No block may differ by more than 64 of 255 levels in any color channel.

Here is how 64 levels was chosen:

- The same panels drawn with a different renderer (Windows GDI instead of cairo) differ by at most 49 levels.
- A real change in the data differs by much more: a point moved by 15 percentage points differs by 188 levels, and a cutoff line moved by 0.1 by 99.
- A panel that positions marks by measuring its own text (for example Figure S1C) can shift as a whole when a font measures differently, and then fails. If that happens, compare the panel with the published one by eye.

## How closely the tables match

Every table is compared with the published one, and every table a panel is drawn from has a published copy. They must match character for character, with two exceptions, each stated in the `reason` column of `expected/<figure>.tsv`:

- **Gene order inside enrichment results.** clusterProfiler lists the genes behind each enriched term in an order that changes between runs of the same code. Those gene lists are compared as sets (Figures 2E and S5). Terms, row order and every number must still match exactly.
- **Row order.** Where the published table's row order carries no meaning, the rows may come in any order (`rows_unordered` in `expected/<figure>.tsv`; Figures 5A and S3A). Every row must still match exactly.

No table is compared within a numeric tolerance.

## Names used inside tables and file names

The tables and deposited files keep the column names, values and file names they were made with, so some use shorter or older terms than the paper. This list gives their meaning.

| Name in tables or files | Meaning |
|---|---|
| `author_label`, `author_label_original`, `author_label_class`, `author_layer`, `AUTHOR`, `primary_label` in file names | the published label: the cell type the original study assigned |
| `gut`, `gut_*`, "Gut EC", "Primary gut EC", "Gut EC or enteroendocrine" | intestinal; primary intestinal EC cells from the intestinal tissue studies |
| "Native intestine" (column `context` of `figure1/endocrine_display_cells.tsv.gz`) | primary intestinal tissue |
| `marker_rule`, `*_marker_rule`, `MARKER_RULE`, `rule_called` | the label-free four-marker EC classifier |
| `HIGH_CONFIDENCE_EC`, `LOW_CONFIDENCE_EC`, `INSUFFICIENT_EC_EVIDENCE` | the four-marker classifier's three classes: passes the score and marker-detection cutoffs and the neighborhood check (classified as EC); passes the cutoffs only; fails the cutoffs |
| `candidate`, `four_marker_candidate` | passes the four-marker score and marker-detection cutoffs, before the neighborhood check |
| `veto`, `vetoed`, `local_q`, `locally_consistent` | the neighborhood check: a candidate EC cell is kept only if other candidates are enriched among its nearest neighbors in expression (q ≤ 0.05) |
| `shipped27` | 27 enteroendocrine and islet hormone and marker genes left out of the neighborhood graph in one of the compared settings (Figure S2 tables) |
| `heldout`, `held_out`, `heldout_studies/` | test studies: studies the classifier was applied to without being fitted on them |
| `TRANSPORTED_NOT_IN_FIT` | cells scored with the fitted classifier that were not used to fit it (the HIO cells) |
| `family`, `study_family`, `family_id` | study |
| `unit`, `biological_unit` | the replicate: a donor, or a sequenced sample where donors are not identified |
| `arm`, "Leave one arm out" | a group of cells analyzed together (one study, or one condition within a study) |
| `in_dish`, "in dish" | in vitro (stem cell-derived islets in culture) |
| `thinned`, `unthinned` | with and without binomial downsampling to equal sequencing depth |
| `SF_*` (for example `SF_YU2021_ENDODERM`) | study codes; `config/study_names.tsv` gives the study name for each |
| `bank` | the regulon network a regulon comes from |
| `panel5_logFC`, `panel5_FDR`, `paneld_*`, `universe_paneld_tested` | values from the Figure 2C donor-level fit, and the genes it tested; `panel5_core_welch` is always FALSE |
| `old_rule_D_show_in_5E`, `show_in_5E` | whether the factor is selected and more than half of its passing networks overlap G2 or G3; reported only, not used to choose what Figure 5B draws |
| `standardisation`, `neighbourhood`, `paediatric` | British spellings kept in parameter files and study keys |
| prefixes such as `A_`, `C_`, `E_`, `c5_`, `f4b_`, `f5_5a_`, `s1_`, `s2_`, `q1_`, `f3_03_`, `fig3_` | groups of files; not figure numbers (see `PANELS.tsv` for the panel each file feeds) |

The gene groups are G1, EC genes low in SC-SICs; G2, shared EC genes; and G3, a separate SC-SIC program (Figure 2D).

## Licenses and citation

Code: MIT license (`LICENSE`). Deposited data: CC BY 4.0. The deposited data derive from the public studies listed in `upstream/SOURCES.tsv`; cite them as listed there. MSigDB gene sets are CC BY 4.0; cite MSigDB and, for any C8 cell-type set, its source study (`upstream/msigdb_c8/`). The two preprint gene lists are not part of the deposit; cite the preprints as listed in `upstream/SOURCES.tsv`.
