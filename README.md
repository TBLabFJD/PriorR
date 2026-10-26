# PriorR v4.0

PriorR is a prioritization program of disease-linked genetic variants developed within the Genetics & Genomics
Department of La Fundación Jiménez Díaz University Hospital. It analyses the output of the FJD pipeline of SNVs and
CNVs: filtering by gene panels, population frequencies, pathogenicity predictors, inheritance patterns,
phenotypes (HPO), regions and variants already found by another test of the same patient.

## What's new in 4.0

- **Any file, the same way**: gene panels, exomes and genomes are opened from the disk (no size limit) with
  load filters that keep only the candidates. Big files are cached: opening them again takes seconds.
- **Filtering profiles**: None (default, nothing filtered), WES and WGS, plus the profiles each user saves
  with a name (filters and displayed columns).
- **Gene panels**: PanelApp (Genomics England and Australia), RetiGene (inherited retinal diseases), Orphanet,
  GLOWgenes disease / class specific genes and the panels saved by the users.
- **Trios**: de novo filter (trio chosen among the samples of the file or a PED file).
- **ACMG**: criteria suggested from the annotation, confirmed by the analyst and saved with the
  classification; **reanalysis** of classified variants whose ClinVar, frequency, PanelApp or ACMG
  suggestion changed.
- **Scores**: filters for every score of the file (REVEL, AlphaMissense, SpliceAI, ...); AlphaMissense
  and the gnomAD gene constraint are added when the file doesn't have them.
- **GLOWgenes**: view of candidate variants ordered by the GLOWgenes rank of their gene (column of the file,
  Orphanet and PanelApp rankings, or rankings computed by the lab with the GLOWgenes Docker image), and the
  gene-disease specificity score (SGDS).
- **Users and groups**: comments and classifications of variants shared with the user's groups.
- **Reproducible analyses**: the same filters always give the same table and prioritization.
- **Fast tables**: processed in the server (100,000 rows in ~3 s), fixed column widths.
- **Reads**: BAM / CRAM reads of the selected variant with igv.js (files stay on the computer), and
  links to IGV desktop, GeneCards, STRING, gnomAD, VarSome, SpliceAI...
- R 4.4 only (no Python), modular code, automated tests. The VCF annotation of PriorR 3 was removed.

## Installation (shared folder, each user runs it on their computer)

1. Put this folder in a shared folder that every computer mounts (e.g. `\\server\genetica\PriorR`).
2. Install **R 4.4 or newer** on each computer (https://cran.r-project.org).
3. Start PriorR with **`PriorR.bat`** (Windows), **`PriorR.command`** (macOS) or **`priorr.sh`** (Linux).
   It opens in the browser. The first time on each operating system it installs the R packages it needs in
   `library/` (a few minutes); the next computers with the same system don't install anything.
4. The first user creates the administrator.

The shared data stays in `shared/` next to the code:

```
shared/accounts.json           users and groups (passwords encrypted)
shared/comments/<user>.jsonl   comments and classifications of each user
shared/profiles/<user>.json    filtering profiles of each user
shared/sessions/<user>.jsonl   saved sessions of each user
shared/bookmarks/              state of the saved sessions
shared/panels/                 gene panels saved by the users
```

Each user only writes their own files, so several computers can use it at the same time over the network.
Group privacy is enforced by PriorR: anyone with access to the shared folder can read its files.

The WGS cache is kept in the local disk of each computer (`PRIORR_CACHE_DIR`).

PriorR can also run as a server for everybody (one computer, users connect with the browser) with Docker:

```bash
docker build -t priorr .
docker run -p 8888:8888 -v priorr:/var/priorr -v /path/to/variant/files:/data priorr
```

In that case serve it behind HTTPS.

## Using it

- **Load variants > Open variant file**: choose the file, tick a profile (None, WES, WGS or one of yours)
  and check the load filters. The profile also sets the filters of the sidebar and the displayed columns.
- **Profiles**: apply a profile, or save the current filters and columns as a profile of yours.
- **Genes and panels**: filter by panels of the library (search by name), include PanelApp amber genes,
  upload a panel and save it in the library, exclude the genes of other panels.
- **Classify** (button in each row): comment and classify a variant, shared with all your groups, one group
  or nobody. **Classification > Show classifications and comments** adds them to the table.
- **Users and groups** (administrators): users, passwords, groups, and the PanelApp update.

Command line (from the `app` folder):

```bash
Rscript Program_files/scripts/priorr_admin.R create-admin <user>      # asks for the password
Rscript Program_files/scripts/priorr_admin.R reset-password <user>
Rscript Program_files/scripts/priorr_admin.R import-v3 <user> [file]  # classifications of PriorR 3
Rscript Program_files/scripts/update_panels.R                         # download PanelApp
Rscript Program_files/scripts/import_retigene.R <TableS1.xlsx>        # RetiGene (Table S1 of the paper)
```

## Configuration

Settings are read from `priorr.cfg` files (`key: value` lines, `#` for comments): first the one of
the computer (written by the first-run setup), then the one of the lab next to `launch.R`
(see `priorr.cfg.example`): `shared_folder`, `data_folders`, `cache_folder`, `memory_limit`,
`threads`, `genome_build`. Environment variables override them:

| Variable | Default | |
|---|---|---|
| `PRIORR_SHARED_DIR` | `shared` | Users, groups, comments, profiles, sessions, saved panels |
| `PRIORR_DIR` | `app/Program_files` | Data shipped with PriorR (panels, HPO, profiles, scripts) |
| `PRIORR_DATA_DIRS` | home, PriorR folder and disks | Folders "Open file" can browse: `Name=path;Name2=path2` |
| `PRIORR_CACHE_DIR` | local user cache | Cache of the big variant files |
| `PRIORR_WGS_MEMORY` | `2GB` | Memory limit while reading big files |
| `PRIORR_WGS_THREADS` | `4` | Threads while reading big files |
| `PRIORR_PORT` | random | Port of the local PriorR |

Built-in profiles are JSON files in `app/Program_files/profiles` (they can be edited or added) and the
links of the variant details are in `app/Program_files/links.tsv`. The in-app **Help**
documents the installation for other labs, the configuration, the profile format and the reads viewer.

## Code

```
launch.R, PriorR.bat, PriorR.command, priorr.sh   start PriorR on a computer
app/app.R                  entry point (app/R/*.R are loaded automatically)
app/R/snv_filter.R         all the SNV filters (pure function: data + parameters -> table)
app/R/cnv_filter.R         CNV filters
app/R/wgs_prefilter.R      reading of big files with load filters (DuckDB, Parquet cache)
app/R/prioritization.R     red / orange prioritization, Orphanet GLOWgenes ranking
app/R/inheritance.R        inheritance patterns and ROH
app/R/profiles.R           filtering profiles
app/R/panels.R             gene panel library and PanelApp download
app/R/store.R              shared data (users, groups, comments, sessions)
app/R/mod_*.R              Shiny modules
tests/testthat             automated tests (synthetic data)
```

Run the tests from the repository root: `Rscript tests/run_tests.R`

## Developers

- Raquel Romero Fernández (main developer of PriorR 1-3)
- Gonzalo Núñez Moreno
- Pablo Mínguez (PriorR 4)

## License

PriorR source code is provided under the [**Creative Commons Attribution-NonCommercial-ShareAlike 4.0 International (CC BY-NC-SA 4.0)**](https://creativecommons.org/licenses/by-nc-sa/4.0/). PriorR includes several third party packages provided under other open source licenses, please check them for additional details. Gene panels from [PanelApp](https://panelapp.genomicsengland.co.uk) (Genomics England), [PanelApp Australia](https://panelapp-aus.org) and [RetiGene](https://retigene.erdc.info) (Am J Hum Genet 2025;112(10):2253-2265, doi:10.1016/j.ajhg.2025.08.017, CC BY 4.0). GLOWgenes rankings of
PanelApp panels, specific genes and SGDS from Benítez Y, Uría-Regojo G, Mínguez P. A global survey of systems biology-based predictions of gene-rare disease associations to enhance new diagnoses. Sci Rep 2026. doi:10.1038/s41598-026-54510-3 (CC BY-NC-ND 4.0, unmodified, see `app/Program_files/GLOWgenes`).
Human Phenotype Ontology (https://hpo.jax.org, `app/Program_files/HPO`, updated with
`Program_files/scripts/update_hpo.R`).
Gene constraint of gnomAD v4.1 (https://gnomad.broadinstitute.org, `app/Program_files/Constraint`).
Pathogenic missense variants of ClinVar (Landrum MJ et al., Nucleic Acids Res 2018; https://www.ncbi.nlm.nih.gov/clinvar,
`app/Program_files/ClinVar`, updated with `Program_files/scripts/update_clinvar_missense.R`).
AlphaMissense (Cheng et al., Science 2023; https://zenodo.org/records/10813168, listed there as
CC BY 4.0, the header of the file says CC BY-NC-SA 4.0): not included, downloaded by the lab with
`Program_files/scripts/get_alphamissense.R`.
PriorR is developed by the [Bioinformatics Unit - Genetics and Genomics Department, IIS-Fundación Jiménez Díaz](https://www.translationalbioinformaticslab.es/).

[![Licencia de Creative Commons](https://i.creativecommons.org/l/by-nc-sa/4.0/88x31.png)](http://creativecommons.org/licenses/by-nc-sa/4.0/)
