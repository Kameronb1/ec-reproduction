# Enterochromaffin cell analysis: reproducing the figure panels

This repository recreates the main and supplementary figure panels of the paper, from Figure 1C on, starting from deposited input data. For each panel it runs four steps:

1. The analysis script computes the panel's numbers.
2. Those tables are compared with the published source data.
3. The plot script draws the panel from those checked tables.
4. The drawn image is compared with the published panel.

Figures 1A and 1B are diagrams. Their tables are recreated and checked, but they are not drawn.

## What you need

- The input data deposit (about 5 GB): Zenodo DOI [10.5281/zenodo.24093050](https://doi.org/10.5281/zenodo.24093050). Unpack it into a folder. One input is downloaded from its original source instead: MSigDB 2026.1 (Zenodo record 18968178). Its URL and checksum are in `inputs/`.
- Python 3.13 and [Snakemake](https://snakemake.readthedocs.io) 9.23.1: `conda env create -f envs/python.yaml`.
- R 4.6.1 with the packages in `envs/r_packages.tsv` (CRAN plus Bioconductor 3.23): `Rscript envs/install_r_packages.R`, then `Rscript envs/check_versions.R`. R must have the cairo graphics device (`capabilities("cairo")` is `TRUE`), and the Arial font must be installed.
- To install the exact version of every package, dependencies included, as in the run that matched every published table and panel: `pip install -r envs/python_lock.txt` for Python, and `renv::restore(lockfile = "envs/renv.lock")` for R (needs the [renv](https://rstudio.github.io/renv/) package).
- Memory: 8 GB is enough. The largest single step uses about 6 GB. Several steps run at once, so on a computer with little free memory use fewer cores (for example `--cores 2`).

## Run it

```bash
python workflow/scripts/common/fetch_url_inputs.py /path/to/unpacked/deposit
snakemake --cores 4 --resources heavy_r=1 --config data_dir=/path/to/unpacked/deposit
```

The first command downloads the input that comes from its original source into the same folder. `heavy_r=1` makes the memory-heavy R steps run one at a time. With `--cores 8`, a full run took 61 minutes on a Windows computer with 12 cores and 32 GB of memory.

Snakemake works through these steps:

1. Check every input file against the checksum in `inputs/<figure>.tsv`. It stops if anything is missing or changed.
2. Run the analysis scripts, which write tables to `results/<figure>/`.
3. Compare each table with the published one in `expected/<figure>/`. A difference stops the run, and the log shows the first differences. Each passing check is recorded in `results/<figure>/checks/`.
4. Draw each panel as a 600 dpi PNG on a white background: `results/<figure>/panels/<panel>.png`. A panel is drawn only after every table its figure reads has passed its check.
5. Compare each drawn panel with the published one in `expected/<figure>/panels/`. A difference stops the run, and the log shows how large it is. Each passing check is recorded in `results/<figure>/checks/panels/<panel>.png.json`.

To build one figure only, run `snakemake --cores 4 --config 'figures=["figure3"]' data_dir=...`. The Figure 2 supplement also reads Figure 2's tables, so building it alone runs Figure 2's analysis too.

## Find the code for a panel

`PANELS.tsv` has one row per panel. Each row follows the panel from its inputs to its image:

- `shows`: what the panel shows.
- `inputs`: the input manifest.
- `rules`: the Snakemake rules file.
- `scripts`: the analysis scripts.
- `tables`: the panel's tables in `results/`: those it is drawn from and the rest of its source data.
- `checked_against_published`: how many of those tables have a published copy to check against ("no table" for a key shared by several panels).
- `plot`: the plot script.
- `image`: the drawn panel.
- `published_image`: the published panel.

For example, to see how Figure 3D is made, find its row. The analysis script it names computes the numbers, and `workflow/scripts/figure3/plot/3D.R` draws them.

`PANELS.tsv` is built from `panels/<figure>.tsv`. After editing those, rebuild it with `python workflow/scripts/common/build_panels_index.py`.

The `<figure>` folders hold these panels:

- `figure1` to `figure5`: the main figures.
- `figure2_supp`: Figures S1–S5.
- `figure3_supp`: Figures S6–S13 and S15.
- `figure4` also holds Figures S14 and S16, and `figure5` Figures S17–S19, drawn by the scripts in their `plot_supp/` folders.

Panel names, plot scripts and panel images carry the paper's figure numbers and panel letters: `S15A.R` draws Figure S15A. Short prefixes in other file and folder names, such as `s1_marker_limma.R` or `figure3_supp/s9/`, only group files and are not figure numbers; `PANELS.tsv` lists the scripts and tables of each panel.

| Folder | Contents |
|---|---|
| `workflow/scripts/<figure>/` | the analysis scripts; `figure2_supp/classifier_lib/` holds the functions of the 20-gene classifier |
| `workflow/scripts/<figure>/plot*/` | one plot script per panel, `helpers.R` with the figure's shared drawing code, and small drawing files (hand-placed label positions, row layouts, the terms chosen for display) |
| `workflow/scripts/common/` | the input check, the table and image comparisons, and the panel-drawing driver |
| `workflow/rules/<figure>.smk` | which inputs each analysis script reads and which tables it writes |
| `inputs/<figure>.tsv` | every input file: checksum, size, what it is, how it was made |
| `panels/<figure>.tsv` | each panel's tables, scripts, plot script, image name and place on the figure page |
| `expected/<figure>/` | the published tables (also the paper's source data) |
| `expected/<figure>/panels/` | the published panels |
| `config/` | run settings, and the study names used in the figures |
| `envs/` | the Python environment, the R package versions, scripts to install and check them, and lock files with the exact version of every package |
| `upstream/` | how the deposited inputs were made from raw public data |

## How a panel is drawn

Each plot script defines `draw_panel()`. It reads only three things: tables in `results/`, the drawing files in its own folder, and `config/study_names.tsv`. The tables are its own figure's; the Figure 2 supplement also reads Figure 2's. The plotting code shows how each panel is drawn: filtering, aggregation, smoothing, error bars and axis transformations.

`workflow/scripts/common/draw_panel.R` opens a device the size of the figure page (7.5 × 10 in at 600 dpi, cairo). It draws the panel at its place on that page and cuts out the panel's box, given in `panels/<figure>.tsv` as x, y, width and height in pixels.

Drawing on the whole page keeps every text and point position the same as in the published figure. A few panels follow how their published pages were made instead. S6–S8, S10, S11 and S14–S16 are drawn on a canvas of their own and then placed on the page. S2A places its cluster numbers on a canvas the size of the panel.

## How closely the images match

The `image_check` setting in `config/config.yaml` picks the comparison:

- **`exact`.** Every pixel must be identical. This holds on the setup the figures were made with: R 4.6.1 on Windows, the cairo PNG device, Arial, and the package versions in `envs/r_packages.tsv`. On that setup all 86 drawn images (81 panels and 5 keys shared by several panels) match exactly.
- **`tolerance`** (the default). Other systems draw text and line edges slightly differently, so pixels shift without the picture changing. Both images are averaged over 24 × 24 pixel blocks (1/25 in). No block may differ by more than 64 of 255 levels in any color channel.

Here is how 64 levels was chosen:

- The same panels drawn with a different renderer (Windows GDI instead of cairo) differ by at most 49 levels.
- A real change in the data differs by much more: a point moved by 15 percentage points differs by 188 levels, and a cutoff line moved by 0.1 by 99.
- A panel that positions marks by measuring its own text (for example Figure 2C) can shift as a whole when a font measures differently, and then fails. If that happens, compare the panel with the published one by eye.

## How closely the tables match

Every table is compared with the published one. Most must match character for character. A few comparisons are looser, each for a stated reason given in the `reason` column of `expected/<figure>.tsv`:

- **Last-digit differences.** Some numbers differ only in the 14th–16th significant digit. This comes from running on a newer R and a different operating system than the original. These tables are compared as numbers, within a relative tolerance of 1e-12 or 1e-10. Text columns must still match exactly, including every EC classification.
- **Gene order inside enrichment results.** clusterProfiler lists the genes behind each enriched term in an order that changes between runs of the same code. Those gene lists are compared as sets. Terms, row order and every number must still match exactly.
- **Row order.** Where the published table's row order carries no meaning, the rows may come in any order (`rows_unordered` in `expected/<figure>.tsv`). Every row must still match exactly.

Some tables a panel is drawn from have no published copy, because the published figure computed them while drawing. In `PANELS.tsv` these panels show fewer checked tables than tables, for example "2 of 3" for Figure 3B. Each such table is computed from tables that are checked, and the image check confirms what is drawn from it.

One panel's tables are not recomputed. Figure S4C's matched-random-gene-set test takes about an hour on one core and was run on a cluster. Its two result tables are part of the deposit and are copied into `results/`. `classifier.R` makes them when given its `matched_random=` and `hio_matched=` outputs.

## Names used inside tables and file names

The tables keep the column names, values and file names they were published with, so some use shorter or older terms
than the paper. This list gives their meaning.

| Name in tables or files | Meaning |
|---|---|
| `author_label`, `author_label_original`, `author_label_class`, `author_layer`, `AUTHOR`, `ORIGINAL_AUTHOR_*`, `AUTHOR_SILENT_*` | the primary study label: the cell type the original study assigned |
| `gut`, `gut_*`, "Gut EC", "Primary gut EC", "Gut tissue" | intestinal; native intestinal EC cells from the intestinal tissue studies |
| `GUTG`, `GUTT`, `SCG`, `SCT` (Figure S6) | genes chosen in intestinal (`GUTG`) or SC-islet (`SCG`) cells; classifier trained on intestinal (`GUTT`) or SC-islet (`SCT`) cells |
| `MARKER_RULE`, `marker_rule_*`, `rule`, `rule_called`, "marker-rule call" | the label-free four-marker EC classifier (Figure 2A) |
| `MARKER_RULE_GATED`, `MARKER_RULE_UNGATED` | the four-marker classifier with and without its neighborhood check |
| `HIGH_CONFIDENCE_EC`, `LOW_CONFIDENCE_EC`, `INSUFFICIENT_EC_EVIDENCE` | the four-marker classifier's three classes: passes the score and marker-detection cutoffs and the neighborhood check (classified as EC); passes the cutoffs only; fails the cutoffs |
| `CALLED`, `ABSTAIN`, `NOT_CALLED` | classified as EC, too little evidence to classify, not classified as EC |
| `call`, `called`, `*_calls` | whether a cell is classified as EC (or, in gene tables, whether a gene passes its test) |
| `stage_one_call`, `gate_pass`, `two_stage_call` | the 20-gene classifier's score threshold, its neighborhood check, and both together |
| `veto`, `*_VETO`, `CANDIDATE_BEFORE_VETO` | the neighborhood check: a candidate EC cell is kept only if other candidates are enriched among its nearest neighbors in expression (q ≤ 0.05) |
| `v05_call` | the 20-gene classifier's result at the unrelaxed score threshold |
| `heldout`, `held_out`, `heldout_studies/` | test studies: studies the classifier was applied to without being fitted on them |
| `family` | study |
| `unit` | the replicate: a donor, or a sequenced sample where donors are not identified |
| `arm`, "Leave one arm out" | a group of cells analyzed together (one study, or one condition within a study) |
| `is_tp` (Figure 2B) | TRUE for Singh 2024, whose replicates are days of one induced line rather than donors |
| `in_dish`, "in dish" | in vitro (stem-cell-derived islets in culture) |
| `thinned`, `unthinned` | with and without binomial downsampling to equal sequencing depth |
| `shipped27` | 27 enteroendocrine and islet hormone and marker genes left out of the neighborhood graph in one of four compared settings (Figure 2 supplement) |
| "(shipped)" | the training studies as used for the paper |
| `TRANSPORTED_NOT_IN_FIT` | cells scored with the fitted classifier that were not used to fit it (the HIO grafts) |
| `SF_*` (for example `SF_YU2021_ENDODERM`) | study codes; `config/study_names.tsv` gives the study name for each |
| `bank` (Figure 5) | the regulon network a regulon comes from |
| `ALL_ELIGIBLE` | every intestinal study that could be scored agrees |
| `panel5_logFC`, `panel5_FDR`, `paneld_*` (Figure 3) | values from the Figure 3D donor-level fit; `panel5_core_welch` is always FALSE |
| `old_rule_D_show_in_5E`, `show_in_5E` | whether the factor is selected and more than half of its passing networks overlap G2 or G3; reported only, not used to choose what Figure 5E draws |
| `standardisation`, `neighbourhood`, `paediatric` | British spellings kept in parameter files and study keys |
| `*_contract.json`, `classifier_lib/vendor/` | the classifier's fixed parameter files and its shared function library |
| prefixes `s1_` to `s9_`, `f5_s4_`, `f5_5a_` | groups of files; not figure numbers (see `PANELS.tsv` for the panel each file feeds) |

## Licenses and citation

Code: MIT license (`LICENSE`). Deposited data: CC BY 4.0. The deposited data derive from the public studies listed in `upstream/SOURCES.tsv`; cite them as listed there. MSigDB gene sets are CC BY 4.0; cite MSigDB and, for any C8 cell-type set, its source study (`upstream/msigdb_c8/`).
