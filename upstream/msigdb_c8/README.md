# MSigDB C8 cell type signature gene sets

The cell-type enrichments of Figure 3G and Figures S8B and S13B test gene lists against MSigDB collection C8, release 2026.1.Hs.

## Source

- Archive: `msigdb.2026.1.zip` from Zenodo record 18968178, `https://zenodo.org/records/18968178/files/msigdb.2026.1.zip?download=1`. MD5 `512ba99c6827141a9d471972b812d4ac`, SHA-256 `84c8ff3270db946ced2a7a4beacf866a3d4ce138af2cbdb1dae76eacb4d4c576`, 44,335,790 bytes. This URL and MD5 are the ones built into msigdbr 26.1.1.
- Reader: R package msigdbr 26.1.1 (CRAN, MIT license for the package code), called as `msigdbr(species = "Homo sapiens", collection = "C8")`. With human source data no cross-species mapping happens; the export script checks this.
- Members used from the archive: `msigdb.2026.1.Hs.C8.rds` (SHA-256 `753e8355…de0d9`) and `msigdb.2026.1.summary.rds` (SHA-256 `f943c4ab…a7f63`). Both are byte-identical inside the zip.
- Acquisition: `workflow/scripts/common/fetch_url_inputs.py` downloads the URL above and refuses a SHA-256 mismatch; `enrich_c8.R` also checks the MD5 before unpacking.

## Plain-text export

`export_c8_members.R` reads the msigdbr cache and writes the C8 gene sets as plain text:

```
Rscript export_c8_members.R <cache_root> <out_dir>
```

`<cache_root>/R/msigdbr/` must hold `msigdb.2026.1.zip` and its two .rds members. Before it writes any file, the script stops unless:

- the three cache files match the SHA-256 values above;
- msigdbr's output and the raw cached data frame hold the same multiset of (set, symbol) pairs;
- msigdbr did not remap symbols or modify the cache.

It writes three files:

| File | SHA-256 | Content |
|---|---|---|
| `c8_members_msigdb_2026.1.Hs.tsv` | `bed98a6e75060e1ead2e9155eb3ffe3695a37bb1e352eb9b42a5dd6cdc79ae2e` | the exact table the enrichment receives (`gs_name`, `gene_symbol`, in msigdbr's row order and multiplicity) plus set and gene identifiers: 157,637 rows, 866 sets, 20,573 symbols |
| `c8_set_citations_msigdb_2026.1.Hs.tsv` | `6631a4a226218543ab8dc8ff8f70647b31c1e6a040545dbb27b3df8fc8ceb9d7` | one row per set: PMID, GEO ID, exact source table, gene counts |
| `export_summary.tsv` | | row, set and symbol counts, and the R and package versions |

The checksums are those of an export under R 4.6.1 on Windows, with msigdbr 26.1.1, data.table 1.18.6.1 and babelgene 22.9. The citations table is in this folder. The members table (about 12 MB) is not; the script above makes it.

175 rows of the members table repeat a (set, symbol) pair because several source genes map to one symbol. They are kept because the enrichment receives them. All 866 sets carry a source PMID.

Every C8 term tested or called in Figure 3G and Figures S8B and S13B is in this export.

## Terms and attribution

MSigDB gene sets from v2022.1 onward are distributed under CC BY 4.0 (https://www.gsea-msigdb.org/gsea/msigdb_license_terms.jsp, read 2026-10-15). Additional terms on that page apply only to KEGG, BioCarta and the retired STKE sets, none of which are in C8. Attribution: cite MSigDB (Liberzon et al. 2015; Subramanian et al. 2005), msigdbr, and for any reported C8 set its source study, using the PMID in `c8_set_citations_msigdb_2026.1.Hs.tsv`.

## Not covered by these checks

- Whether an msigdbr cache made on another machine gives the same files. Run on such a cache, `export_c8_members.R` first compares its three files with the SHA-256 values above and stops if any differs.
- Equality with Broad's `c8.all.v2026.1.Hs.symbols.gmt`. Not used and not compared.
- The enrichment itself. These checks cover its gene-set input only.
