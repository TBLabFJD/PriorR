# PriorR v4.0 help

PriorR is developed by the [Bioinformatics Unit - Genetics and Genomics Department, IIS-Fundación
Jiménez Díaz](https://www.translationalbioinformaticslab.es/). Code: [github.com/TBLabFJD/PriorR](https://github.com/TBLabFJD/PriorR).

PriorR prioritizes genetic variants from the files of the FJD pipeline (SNVs and CNVs): filters by
gene panels, population frequencies, pathogenicity predictors, inheritance, phenotypes (HPO),
regions, and variants already found by another test of the same patient.

**Contents:** [How it works](#how-it-works) · [Installation](#installation) ·
[Configuration](#configuration) · [Opening a file](#opening) · [Filtering profiles](#profiles) ·
[The profile file](#profile-file) · [Gene panels](#panels) · [Filters](#filters) ·
[Inheritance, family and ROH](#inheritance) · [Prioritization](#prioritization) ·
[GLOWgenes ranking](#glowgenes) · [Variant information](#variant-information) · [Reads (BAM / CRAM)](#reads) ·
[ACMG criteria](#acmg) · [Comments and classifications](#comments) · [Reanalysis](#reanalysis) ·
[Analysis rounds](#rounds) · [Saved sessions](#sessions) ·
[Users and groups](#users) · [CNVs](#cnv) · [Files the lab can adapt](#files) ·
[Troubleshooting](#troubleshooting) · [Citations](#citations)

The complete user guide with screenshots is in the **PriorR 4.0 user guide (PDF)** (link in the box on the right).

<h2 id="how-it-works">How it works</h2>

- **Everybody runs PriorR on their own computer.** It opens in the browser, but it is not a web
  site: the program runs on the computer itself.
- **The variant files (genomes, exomes, panels) and the BAM / CRAM files stay where they are**,
  usually on each person's computer. PriorR reads them from there, without copying or uploading them.
- **The lab shares a folder** (for example on the network drive) where PriorR keeps what is
  common: users, groups, comments and classifications, profiles, saved sessions and the lab's own
  panels. Every computer of the lab must use the same folder.

<h2 id="installation">Installation</h2>

**With the installers** (`PriorR-4.0-setup.exe` for Windows, `PriorR-4.0.pkg` for macOS): wizards that
install PriorR for the user (no administrator needed) and add a PriorR icon; the Windows one also asks
for the shared folder of the lab. R must be installed first: the installers say so when it is missing.
They are not signed: Windows shows "Windows protected your PC" (*More info > Run anyway*) and macOS needs
right click > *Open* the first time.

**Without installers, on each computer:**

1. Install **R 4.4 or newer** from https://cran.r-project.org (Windows, macOS or Linux).
2. Have a copy of the PriorR folder (it can be on the shared drive, so that everybody uses the same
   version, or on each computer).
3. Open PriorR with **PriorR.bat** (Windows), **PriorR.command** (macOS) or **priorr.sh** (Linux).
   The first time on each operating system it installs the R packages it needs in the `library`
   folder of PriorR (a few minutes); afterwards it starts in seconds. Closing the console window
   (the black one) stops PriorR.

**The first time PriorR opens on a computer** it asks where the shared folder of the lab is
(unless the lab configured it, see below). The first person of the lab also creates the
**administrator user**, who adds the others.

**To set up PriorR for another lab**, the easiest way is:

1. Put the PriorR folder on the shared drive of the lab.
2. Create in it the file `priorr.cfg` (a copy of `priorr.cfg.example`) with the shared data
   folder. Nobody has to configure anything: they only install R and open PriorR.
3. Open PriorR, create the administrator and add the users and groups.

<h2 id="configuration">Configuration</h2>

The configuration is in `priorr.cfg` text files with `key: value` lines (lines starting with `#`
are comments). PriorR reads, in this order:

1. The `priorr.cfg` **of this computer** (PriorR writes it the first time; the box on the right
   shows where it is).
2. The `priorr.cfg` **of the lab**, next to `launch.R` in the PriorR folder.
3. The defaults.

| Key | What it is for | Default |
|---|---|---|
| `shared_folder` | Shared folder of the lab (users, comments, profiles, sessions, own panels) | folder `shared` inside PriorR |
| `data_folders` | Folders shown first by "Open file": `Genomes=/data/genomes; Exomes=D:\exomes` | none |
| `cache_folder` | Local folder for the cache of big files | the user's cache folder |
| `memory_limit` | Maximum memory to read big files | `2GB` |
| `threads` | Processors to read big files | `4` |
| `max_rows` | Safety limit of rows loaded from a file (above it nothing is loaded) | `2000000` |
| `genome_build` | Genome build of the files, for the reads viewer and the external links (`GRCh38` or `GRCh37`) | `GRCh38` |

Examples of `shared_folder`: `\\server\genetics\priorr_data` or `Z:\priorr_data` (Windows),
`/mnt/genetics/priorr_data` (Linux), `/Volumes/genetics/priorr_data` (macOS). In the lab `priorr.cfg` it
can also be relative to the PriorR folder (`shared_folder: shared`): then copying the PriorR folder copies
the program and the data of the lab (users, comments, classifications, rankings, AlphaMissense).

Content of the shared folder:

```
accounts.json          users and groups (encrypted passwords)
comments/<user>        comments and classifications of each user
profiles/<user>        filtering profiles of each user
sessions/<user>        saved sessions of each user
bookmarks/             state of the saved sessions
panels/                panels saved by the users
```

Each user only writes their own files, so several computers can use PriorR at the same time.
**Privacy:** PriorR only shows each user the comments they are allowed to see, but anybody with
access to the shared folder can read its files directly.

<h2 id="opening">Opening a variant file</h2>

**Load variants > Open variant file** opens a dialog with four steps:

1. **File:** choose the file of the pipeline (tab-separated text: `.tsv`, `.txt`, also compressed `.gz`)
   from the disk or the network. The dialog lists every file of the folder; a file that is not a
   variant table gives a message when it is loaded. There is no size limit: a 7 GB genome opens like a panel.
2. **Filtering profile:** tick a profile (None, WES, WGS or one of yours). The profile fills the
   load filters and, once loaded, also the sidebar filters and the columns.
3. **Genes and panels:** load only the typed genes, the genes of some panels or of a gene file,
   and leave out the genes of other panels.
4. **Load filters:** only the variants passing them are loaded (empty = no filter). The variants of
   "Always keep" are always loaded: ClinVar pathogenic and high SpliceAI, even outside the regions.

The **first time** a big file is opened PriorR prepares it (1-2 minutes for a genome) and keeps an
optimized copy in the cache of the computer; the next times, with other filters, it takes
seconds. To protect the memory of the computer there is a safety limit of rows (`max_rows` in
`priorr.cfg`, 2,000,000 by default): when more rows pass the load filters, nothing is loaded and PriorR
says how many there are (it never loads only part of them).

<h2 id="profiles">Filtering profiles</h2>

A profile keeps **a whole type of analysis**: the load filters, the sidebar filters and the
displayed columns.

- **None:** no filters (default).
- **WES:** exome or panel: rare (frequency ≤ 0.01) exonic and splicing variants, canonical
  transcript, CADD ≥ 15, no ClinVar benign; always ClinVar pathogenic and SpliceAI ≥ 0.5.
- **WGS:** genome: like WES without the CADD and canonical transcript filters, with up to 200,000 rows.

Each user can **save profiles** with a name (in the open file dialog or in **Profiles** of the
sidebar), **download** them as a `.json` file to give them to a colleague, and **upload** the
profile a colleague gives them. The built-in profiles are files in `app/Program_files/profiles`
and the administrator can change them or add others.

<h2 id="profile-file">The profile file (.json)</h2>

```json
{
  "name": "Retina WES",
  "description": "Exomes of retinal dystrophies",
  "load": {
    "regions": ["EXONIC", "SPLICING"],
    "canonical_only": true,
    "max_af": {"gnomADg_AF_popmax": 0.01, "gnomADg_AF": 0.01, "gnomADe_AF": 0.01, "FJD_MAF_AF": 0.01},
    "min_cadd": 15,
    "exclude_clinvar_benign": true,
    "rescue_clinvar_pathogenic": true,
    "rescue_spliceai_min": 0.5,
    "genes": [],
    "panels": ["retigene:all"],
    "amber": false,
    "excluded_panels": []
  },
  "filters": {
    "region": ["EXONIC", "SPLICING"],
    "canonical": true,
    "max_gnomadg_af": 0.01, "max_gnomadg_popmax": 0.01, "max_gnomade_af": 0.01, "max_fjd_af": 0.01,
    "min_cadd": 15,
    "exclude_benign": true,
    "inheritance": "",
    "panels": [], "amber": false, "excluded_panels": [],
    "prioritize": true, "color_priority": true
  },
  "columns": ["SYMBOL", "HGVSc", "HGVSp", "Consequence", "CLNSIG", "CADD_PHRED", "gnomADg_AF_popmax"]
}
```

**load** (filters when opening the file):

| Field | Meaning |
|---|---|
| `regions` | Genomic regions (`EXONIC`, `SPLICING`, `INTRONIC`, `UPSTREAM`, `DOWNSTREAM`, `3UTR`, `5UTR`, `ncRNA`, `regulatory`); `[]` = all |
| `canonical_only` | Only rows of the canonical transcript |
| `max_af` | Maximum frequency per column; variants without frequency are kept; `null` = no filter |
| `min_cadd` | Minimum CADD PHRED; variants without CADD are kept |
| `exclude_clinvar_benign` | Leave out ClinVar benign / likely benign |
| `rescue_clinvar_pathogenic` | Always load ClinVar pathogenic variants |
| `rescue_spliceai_min` | Always load variants with SpliceAI DS max ≥ this value (and rare) |
| `genes`, `panels`, `amber`, `red` | Load only these genes / the genes of these panels (with PanelApp amber genes if `amber`, amber and red if `red`) |
| `excluded_panels` | Don't load the genes of these panels |

**filters** (sidebar; what is not there doesn't filter): `type`, `region`, `consequence`
(lists), `canonical`, `exclude_predicted` (true / false), `max_gnomadg_af`, `max_gnomadg_popmax`,
`max_gnomade_af`, `max_fjd_af` (0 to 1), `min_cadd` (0 to 75), `scores` (score filters, e.g.
`{"REVELScore": 0.5, "AlphaMissense": 0.564}`), `scores_mode` (`all` or `any`), `exclude_benign`, `inheritance`
(`De novo`, `Dominant`, `Recessive-Homozygous`, `Recessive-Compound heterozygous`, `X-linked` or empty),
`panels`, `amber`, `excluded_panels`, `phenotypes` (HPO identifiers), `glow` (max GLOWgenes rank; -1 = no limit), `color_priority`,
`prioritize`.

**columns:** displayed columns, with the names of the file (e.g. `SYMBOL`, shown as GENE);
`null` = the default ones. CHROM, POS, REF and ALT are always shown.

Panels are written as `source:identifier`: `panelapp_genomics_england:307`,
`panelapp_australia:3149`, `retigene:rp`, `orphanet:retinal_dystrophy_panel`,
`custom:Name of the panel`. The identifier of each panel can be seen in a saved profile that uses it.

<h2 id="panels">Gene panels</h2>

- **PanelApp** (Genomics England and Australia): green (diagnostic) genes; "Include amber genes" adds the
  amber ones and "Include amber and red genes" also the red ones (low evidence). Red genes are only in the
  PanelApp copies downloaded with PriorR 4.0 (if the box says the copy has no red genes, download the
  latest PanelApp). PriorR comes with a copy for the whole lab; in *Genes and
  panels*, "Download the latest PanelApp (only for me)" downloads the current version only for you
  (in `panelapp/<user>` of the shared folder; needs internet, 10-15 minutes): the other users keep
  the lab copy. The lab copy is updated by running `Rscript Program_files/scripts/update_panels.R`
  from `app/`.
- **RetiGene** (gene atlas of inherited retinal diseases, AJHG 2025): a panel with all the genes
  (retained in green, candidates in amber) and one per clinical group (RP, LCA, MD, CD/CRD, CSNB,
  Usher, BBS, Joubert...).
- **Orphanet:** panels by disease groups (their GLOWgenes rankings are in *GLOWgenes ranking*).
- **Disease-specific genes (Benítez et al., 2026):** their own list in *Genes and panels*: genes
  predicted by GLOWgenes as specific of one disease (*disease-specific genes*) or of a class of
  diseases (*class-specific genes*). They are added to the genes of the other panels.
- **Own panels:** upload a gene file (one gene per line) and save it with a name; it is available
  for the whole lab.

**Phenotypes / HPO:** choose one or several HPO terms (by name or `HP:` identifier), or paste a list
of identifiers (`HP:0000510, HP:0000556`). The genes of
each term and of all its subterms (Human Phenotype Ontology, `genes_to_phenotype`) work like one more
gene panel: they are added to the genes of the panels chosen in *Genes and panels*. They are saved in
the profiles (`phenotypes`). The administrator updates HPO with
`Rscript Program_files/scripts/update_hpo.R` (from `app/`, needs internet).

The **SGDS** column (gene-disease specificity score, Benítez et al. 2026) is added to the variants:
close to 1 when the GLOWgenes predictions of the gene concentrate in few diseases (specific gene),
close to 0 when the gene is predicted for many diseases (e.g. highly connected genes).

<h2 id="filters">Sidebar filters</h2>

They apply to the loaded variants and can be changed at any time. Hiding columns doesn't change
the result. Variants without gene, consequence or region are kept (ticking "NA" in the lists).
Frequencies or scores that are empty or not a number count as missing: a missing frequency or CADD
never removes a variant.

**Pathogenicity** has a filter for each score the file has: CADD, REVEL, AlphaMissense, MutScore,
SpliceAI (max delta score), dbscSNV ADA / RF, number of pathogenic predictors, DOMINO, phyloP and the
gene LOEUF of gnomAD (maximum). Only the scores of the loaded file are shown. "All must pass": a
variant without a score passes that filter (e.g. REVEL only exists for missense variants); "at least
one must pass": the variant needs one of the chosen scores passing (e.g. REVEL ≥ 0.5 *or* SpliceAI ≥
0.2).

**Annotations added by PriorR** (only when the file doesn't have them):

- **AlphaMissense** (Cheng et al., Science 2023) score and class of the missense SNVs
  (`AlphaMissense`, `AlphaMissense_class`; likely pathogenic > 0.564, likely benign < 0.34). PriorR
  reads them from `annotations/AlphaMissense_hg38.parquet` of the shared folder (or the
  `alphamissense_file` setting), made once by the administrator with
  `Rscript Program_files/scripts/get_alphamissense.R` (from `app/`; downloads ~640 MB from Zenodo).
- **Gene constraint of gnomAD v4.1:** `gnomAD_LOEUF` (< 0.6: intolerant to loss of function),
  `gnomAD_pLI` and `gnomAD_mis_z` (`Program_files/Constraint`, updated with
  `Program_files/scripts/update_constraint.R`).
**Regions (BED)** uses BED coordinates (0-based start). Column names are shown without `_` and
`SYMBOL` as **GENE**; the downloaded files keep the original names.

<h2 id="inheritance">Inheritance, family and ROH</h2>

Without a family the genotype of the sample is used. In *Family / ROH* choose the family:

- **trio (samples of the file):** the proband (affected), their sex, the father and the mother among
  the samples of the file (columns `<sample>_GT`).
- **PED file:** any family. Format: one line per person, separated by tabs or spaces:

```
family  individual  father  mother  sex  phenotype
F1      P001        P002    P003    1    2
F1      P002        0       0       1    1
F1      P003        0       0       2    1
```

Sex 1 = male, 2 = female; phenotype 2 = affected, 1 = unaffected; 0 = unknown. The individual must
match the name of the columns of the file (`P001_GT`).

- **De novo:** needs the proband and both parents. The proband is heterozygous (VAF 0.2-0.8) and both
  parents are homozygous reference with depth ≥ 10 and at most 1 read with the variant (when the file
  has `_DP`, `_AD`, `_VAF`). Chromosome X in a male proband: only the mother is checked (he is
  hemizygous); chromosome Y: only the father. A parent without a call (`./.`) never gives a de novo.
  The thresholds are `denovo_settings` in `app/R/inheritance.R`.
- **Dominant:** affected heterozygous, unaffected homozygous reference, dominant gene in OMIM.
- **Recessive-Homozygous:** affected homozygous alternative, parents carriers, unaffected not homozygous.
- **Recessive-Compound heterozygous:** at least two different heterozygous variants in the same
  recessive gene; with both parents, one inherited from each.
- **X-linked:** chromosome X, affected males hemizygous, affected females carriers, unaffected males without the variant.

**ROH:** colors the variants in runs of homozygosity (`_ROH` columns) of the proband or of the
affected individuals of the PED.

<h2 id="prioritization">Prioritization</h2>

Above the table, the **PriorR priority** button shows the variants of the sidebar filters with the
red rows at the top, then the orange ones, then the rest. In the **Sidebar filters** view, "Color red /
orange" colors those rows without changing the order:

**Red** (any of these; never when ClinVar says benign, likely benign or uncertain significance):

- ClinVar pathogenic or likely pathogenic.
- Rare (gnomAD genomes popmax < 0.001 or absent), more than 60% of the predictors pathogenic and CADD > 25
  and / or MutScore > 0.8.
- Rare start / stop lost, stop gained, frameshift or in-frame indel with more than 75% pathogenic
  predictors (when there is no CADD or MutScore).
- Rare splicing or intronic variant with dbscSNV ADA or RF > 0.8.

**Orange** (any of these):

- Low frequency (popmax 0.001-0.01), more than 60% of the predictors pathogenic and CADD > 25 or
  MutScore > 0.8.
- Low frequency splicing or intronic variant with more than 60% pathogenic predictors and ADA or RF > 0.8.
- Other pathogenic ClinVar annotations (e.g. low penetrance, risk factor), not conflicting.
- Rare ClinVar uncertain significance with more than 60% pathogenic predictors and CADD > 25 or
  MutScore > 0.8.

A variant matching both rules is red. The same filters always give the same table. In the PriorR
priority view, the box above the table says how many red and orange variants there are and repeats
these rules.

<h2 id="glowgenes">GLOWgenes ranking</h2>

The button **GLOWgenes ranking** above the table shows only candidate variants ordered by the
GLOWgenes rank of their gene: rank 0 are the genes of the disease panel and then the closest genes
(1, 2, ...). Inside each gene, PriorR red / orange variants first and then by CADD. One row per
variant (the canonical transcript when there is one).

In *GLOWgenes ranking* of the sidebar, "Rank the genes with":

- **The GLOWgenes column of the file.**
- **Precomputed GLOWgenes rankings** (also used when the file has no `GLOWgenes` column), one or
  several (the best rank of each gene):
  - **Orphanet:** rankings of the Orphanet panels of PriorR.
  - **PanelApp (GLOWgenes, Sci Rep 2026):** rankings of 209 PanelApp panels from Benítez et al. 2026
    (grouped by disease class).
  - **Computed by the lab:** rankings made with *run GLOWgenes* (below), shared with the whole lab.
- **A new ranking (run GLOWgenes):** choose the genes (a panel of the library or a file, at least 10
  genes) and a name, and PriorR runs GLOWgenes with Docker (image `tblabfjd/glowgenes`) in the background; the result
  is saved in the shared folder (`glowgenes/<name>.tsv`). It needs Docker on the computer; the first
  run downloads the image (~1 GB, 4 GB on disk) and a run can take a long time. When it finishes, the
  ranking is added to the precomputed ones and selected.
- **Only genes up to a GLOWgenes rank:** keeps only the genes with a rank up to the chosen one (with
  the chosen rank source), in the GLOWgenes view and in the sidebar view.
- **Filters** (instead of the sidebar ones): rare variants (frequency ≤ 0.01), not ClinVar benign,
  and at least one evidence of pathogenicity: ClinVar pathogenic, high impact (stop, frameshift,
  canonical splicing, start lost), CADD ≥ 20, REVEL ≥ 0.5, SpliceAI ≥ 0.2 or PriorR red / orange.
  They are in `app/Program_files/glowgenes.json` and the lab can change them.
- **Hide panel genes** (next to the buttons, in the GLOWgenes view): hides the genes of the disease
  panel (rank 0), so only the new candidate genes are shown.
- **Sidebar filters** goes back to the usual view.

<h2 id="variant-information">Variant information</h2>

Selecting a row shows the impact of the variant (the `IMPACT` column of VEP when the file has it;
otherwise PriorR derives it from the consequence: HIGH = stop, frameshift, canonical splice site,
start / stop lost; MODERATE = missense, in-frame, splice donor region; LOW = synonymous, splice
region; MODIFIER = the rest: intronic, UTR, up / downstream...), consequence and HGVS, **genotypes** of each sample (GT, DP, AD, VAF and
what each variant caller says), **read depth** (the box is the depth of all the variants of the
file and the red dot the depth of this variant: if it is far to the left, the variant has low
coverage), frequencies, predictors and scores (with AlphaMissense and the gene constraint when available),
OMIM, the [ACMG criteria](#acmg) suggested by the annotation, comments and **links**: to the gene (GeneCards, STRING,
OMIM, ClinVar, gnomAD, PanelApp, DECIPHER, GTEx, UniProt, Ensembl, PubMed) and to the variant
(gnomAD, VarSome, ClinVar and dbSNP when it has an rsID, SpliceAI Lookup, UCSC). The list of
links is in `app/Program_files/links.tsv`. Long values of the table are shown in full on hover
or click.

<h2 id="reads">Reads (BAM / CRAM)</h2>

The **Reads (IGV)** box shows the reads of the selected variant with igv.js (the web version of
IGV): reference, genes, coverage and reads, with zoom and scroll as in IGV.

- **Open BAM / CRAM files:** choose the files of the samples (several at once; the dialog lists the
  `.bam` / `.cram` files, also in capitals). The index must be next to each file: `X.bam.bai`, `X.bai`,
  `X.bam.csi`, `X.cram.crai` or `X.crai` (also in capitals). A file is assigned to the
  sample whose name it contains (`run1_P001.bam` > `P001`).
- The files are read from this computer; PriorR only gives igv.js the pieces it needs. The
  reference genome and the genes come from the IGV servers (internet). CRAM files need the same
  reference genome they were made with (GRCh38 / hg38 by default).
- **Open in IGV desktop:** loads the files in IGV desktop and goes to the variant. IGV must be open
  with the port enabled (*View > Preferences > Advanced > Enable port*, 60151).
- Saved sessions remember the files.

<h2 id="acmg">ACMG criteria</h2>

The **ACMG (suggested)** card of *Variant information* shows the criteria that the annotation
supports (ACMG/AMP 2015 with the ClinGen SVI recommendations), their points and the class
(Tavtigian et al. 2020: very strong 8, strong 4, moderate 2, supporting 1, benign negative;
≥ 10 pathogenic, 6-9 likely pathogenic, 0-5 uncertain, -1 to -6 likely benign, ≤ -7 benign):

- **PVS1:** null variant (stop, frameshift, canonical splice site, start lost). Check that loss of
  function is the mechanism of the disease, the last exon / NMD and the transcripts.
- **PM2_Supporting:** absent or popmax ≤ 0.0001 in gnomAD. **BA1 / BS1:** popmax > 0.05 / > 0.01
  (BS1 depends on the disease).
- **PM4:** in-frame indel or stop lost.
- **PP3 / BP4:** REVEL for missense (Pejaver et al. 2022: PP3 ≥ 0.644 / 0.773 / 0.932, BP4 ≤ 0.290 /
  0.183 / 0.052), CADD when there is no REVEL; SpliceAI ≥ 0.2 (PP3) or ≤ 0.1 (BP4) for splicing.
- **BP7:** synonymous with SpliceAI ≤ 0.1.
- **PS1 / PM5:** the same amino acid change (by another nucleotide change) / another missense change at
  the same residue is pathogenic or likely pathogenic in ClinVar (at least one star), in the same
  transcript. From `Program_files/ClinVar` (`Rscript Program_files/scripts/update_clinvar_missense.R`).
  PS1 needs checking that the variant doesn't change splicing.
- **PM6:** de novo in the trio (PS2 if the parentage is confirmed).

They are only suggestions from the annotation: segregation, functional studies, phenotype, *in
trans* (PM3)... are added by the analyst when classifying.

<h2 id="comments">Comments and classifications</h2>

The **Classify** button of each row classifies the variant (Benign ... Pathogenic) and comments
on it. The suggested **ACMG criteria** come already chosen: remove, change or add criteria and the
points and class are recalculated; the criteria are saved with the classification. Each comment is shared with **all your groups**, **one group** or **nobody**. Each user
sees their own comments and the ones shared with their groups. **Classification > Show
classifications and comments** adds to the table the latest visible classification and the number
of comments.

<h2 id="reanalysis">Reanalysis</h2>

When a variant is classified PriorR saves what the annotation said: ClinVar, gnomAD popmax, the
green PanelApp panels of the gene and the class suggested by ACMG. When a file is opened later
(for example annotated again with a newer ClinVar), a blue box above the table lists the
classified variants of the file that changed since their classification: new ClinVar
classification, frequency crossing 0.0001 / 0.01 / 0.05, the gene in new green PanelApp panels
(of your PanelApp copy) or a different ACMG suggestion. **Classification > Only classified
variants with changes since their classification** shows only those variants.

<h2 id="rounds">Analysis rounds</h2>

To analyze a genome in steps (first a panel, then the phenotype, then everything rare...) without
seeing again what was already reviewed: in **Analysis rounds**, *Save a round with the current
panels, frequency and CADD*. A round keeps the variants of the current gene panels / HPO terms /
genes, maximum frequency and minimum CADD; the other filters (consequence, region, inheritance,
scores...) are not used. Its default name says what it holds, e.g. `P001 · RetiGene - RP (v2025) · AF ≤
0.01 · CADD ≥ 20` (or type another name). With the next filters, tick **Hide the variants of previous
rounds**: only the variants not seen in any round are shown (the line above the table says how
many are hidden). It is optional: untick it to see everything again. The rounds are saved per user
and file in the shared folder (`rounds/`), so they are there when the file is opened again; each
round can be deleted.

<h2 id="sessions">Saved sessions</h2>

**Saved sessions** keeps the current analysis (file, profile, filters, columns, BAM files).
Opening the session loads the file again with the same filters, also from another computer, as
long as the files are reachable at the same path.

<h2 id="users">Users and groups</h2>

The **Users and groups** tab (administrators only) creates users, changes passwords, deactivates
users, makes administrators, creates groups and assigns members. The members of a group see the
comments the others share with their groups.

**Remember me on this computer:** PriorR asks for the password only the first time on each
computer. *Sign out* makes the computer forget the user. When an administrator changes the
password of a user, every computer forgets that user.

Only the first user is created on the sign in screen (the first time a shared folder is used); the
next ones are created by an administrator in **Users and groups**.

From the command line, in the `app` folder of PriorR (it asks for the password twice):

```
Rscript Program_files/scripts/priorr_admin.R reset-password <user>   # new password (e.g. the administrator forgot it)
Rscript Program_files/scripts/priorr_admin.R create-admin <user>     # a new administrator
```

When PriorR uses its own `library/` add `R_LIBS=../library/<platform>` before `Rscript` (e.g.
`R_LIBS=../library/x86_64-pc-linux-gnu-R4.6`).

<h2 id="cnv">CNVs</h2>

The **CNV** tab opens the structural variant files of AnnotSV (`.tsv`) with **Open CNV file**.

- **Annotation mode:** `full` shows one row per CNV (default) and `split` one row per gene of each CNV.
- **Filters:** samples, type (DEL, DUP...), location, minimum ACMG class, minimum number of programs
  (only for files with `N_PROGRAMS`), genes and panels (the same panel library as the SNVs, with the
  disease-specific genes of Benítez et al. 2026), **HPO terms** (their genes work like a panel) and
  maximum frequency of benign gains / losses. By default nothing is filtered. A CNV passes a gene
  filter when any of its genes is in the list.
- **Order by ACMG class and AnnotSV score** puts the most pathogenic CNVs first.
- **GLOWgenes ranking:** choose precomputed rankings (Orphanet, PanelApp, the lab's): the
  `GLOW_rank` column is the best rank of the genes of each CNV; order by it or keep only the CNVs
  with a gene up to a rank.
- **Reads (IGV):** like the SNVs (open the BAM / CRAM files): a small CNV (≤ 20 kb) is shown whole,
  a big one by its two breakpoints (±300 bp), where the reads that support it are.
- The genotype of each sample (`0/1`, `1/1`) is taken from its column of the file.
- **CNV information:** selecting a row shows the position, size, class and AnnotSV criteria, the
  genes, the known pathogenic CNVs of the region, OMIM / GenCC, the constraint of the genes and
  links to UCSC, DECIPHER and gnomAD-SV.

<h2 id="files">Files the lab can adapt</h2>

| File | What for |
|---|---|
| `priorr.cfg` (next to `launch.R`) | Configuration of the lab (see [Configuration](#configuration)); commented example in `priorr.cfg.example` |
| `app/Program_files/profiles/*.json` | Built-in profiles (None, WES, WGS); they can be edited or added with the format of [The profile file](#profile-file). They are listed in the order of their file names |
| `app/Program_files/links.tsv` | Links of the variant and CNV information: type (`gene`, `variant` or `cnv`), name, URL with `{gene}`, `{chrom}`, `{pos}`, `{ref}`, `{alt}`, `{rsid}`, `{ucsc_build}`... and description |
| `app/Program_files/glowgenes.json` | Filters of the GLOWgenes ranking view: `max_af`, `exclude_clinvar_benign`, `evidence` (`clinvar_pathogenic`, `high_impact`, `min_cadd`, `min_revel`, `min_spliceai`, `priorr_priority`; `null` = not used), `one_row_per_variant`, `max_rank` |
| `app/Program_files/Panels/` | PanelApp and RetiGene panels (updated with the scripts, not by hand) |
| `app/Program_files/HPO/` | Human Phenotype Ontology (`update_hpo.R`) |
| `app/Program_files/Constraint/` | Gene constraint of gnomAD v4.1 (`update_constraint.R`) |
| `app/Program_files/ClinVar/` | Pathogenic missense variants of ClinVar for PS1 / PM5 (`update_clinvar_missense.R`) |
| `annotations/` of the shared folder | AlphaMissense (`get_alphamissense.R`) |
| `app/Program_files/examples/` | Synthetic demo data (a trio and CNVs, `make_demo_data.R`) to try PriorR |
| `app/Program_files/GLOWgenes/` | GLOWgenes rankings of PanelApp panels, specific genes and SGDS (Sci Rep 2026, unmodified) |
| `app/Program_files/Paneles_Orphanet/` | Orphanet panels and GLOWgenes rankings |

The changes in these files are seen when PriorR is opened again.

<h2 id="troubleshooting">Troubleshooting</h2>

- **"Maximum upload size exceeded":** don't upload big files; open them with *Open variant file*.
- **The first load of a genome is slow:** it is normal (it prepares the cache); the next ones are fast.
- **"N rows pass the load filters, more than the safety limit":** above `max_rows` (`priorr.cfg`,
  2,000,000 by default) nothing is loaded (PriorR never loads only part of the variants that pass the
  filters). Use stricter load filters (WGS profile, genes, panels) or raise `max_rows`.
- **I don't see the comments of a colleague:** check that you are in the same group and that the
  comment is not private.
- **Each computer sees different users:** they are not using the same shared folder (see the
  configuration box of this page).
- **A BAM / CRAM doesn't appear in the dialog:** check that the file is really in the folder (sometimes
  only the `.bai` or the `.md5` were copied).
- **"Permission denied" when opening a network folder:** PriorR must run as the user who can read it.
  Folders mounted with sshfs (FUSE) only let in the processes of the user who mounted them with their own
  main group: start PriorR from that user's session.
- **The reads don't show:** check that the index is next to the BAM / CRAM, that there is internet
  (for the reference genome) and that `genome_build` matches the files.
- **PriorR doesn't start:** check that R is installed; the first time it needs internet to install
  the packages.

<h2 id="citations">Citations</h2>

- GLOWgenes: de la Fuente L et al. Prioritization of new candidate genes for rare genetic diseases by a
  disease-aware evaluation of heterogeneous molecular networks. Int J Mol Sci 2023;24(2):1661.
- GLOWgenes rankings of PanelApp panels, disease-specific genes and SGDS: Benítez Y, Uría-Regojo G, Mínguez P. A global survey of systems biology-based predictions of gene-rare disease associations to enhance new diagnoses. Sci Rep 2026. doi:10.1038/s41598-026-54510-3.
- RetiGene: Am J Hum Genet 2025;112(10):2253-2265. AlphaMissense: Cheng J et al., Science 2023.
  ClinVar: Landrum MJ et al., Nucleic Acids Res 2018. gnomAD v4.1. Human Phenotype Ontology.
- ACMG/AMP: Richards S et al., Genet Med 2015; Tavtigian SV et al., Hum Mutat 2020; Pejaver V et al., Am J
  Hum Genet 2022; Walker LC et al., Am J Hum Genet 2023.
