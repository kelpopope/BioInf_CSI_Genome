# CSI Genome: The Contamination Detective

**Project 07 · Introduction to Bioinformatics · Astana IT University**  
**Abylaikhan Torekhan & Kuanysh Bakizhan · BDA-2405**

All analysis code is in **[CSI_Genome_Project07.ipynb](CSI_Genome_Project07.ipynb)**.
Its eight tasks follow the handbook in order. No separate analysis scripts are required.

The full run screens **3,006 real GenBank Mollicutes assemblies**.
**246 assemblies (8.18%) are flagged:** 243 have known adapter motifs and three have strong foreign-DNA matches.
NCBI BLAST supports E. coli in one case and mouse mitochondrial DNA in two.
These are screening flags, not a measurement of all contamination in GenBank.

## Read the project

- [Notebook with the code and executed full-run outputs](CSI_Genome_Project07.ipynb).
- [Notebook as HTML](notebook.html), for reading without Jupyter.
- [Report PDF](report/CSI_Genome_Report.pdf) and [editable Word report](report/CSI_Genome_Report.docx).
- [Ten curator candidates](results/curator_top10.csv) and [case review](results/review_decisions.csv).

## Open and run the notebook

Use **Python 3.12**. A C compiler is needed for mappy (macOS Command Line Tools or Linux build-essential).
On macOS/Linux, run from the repository folder:

```sh
python3.12 -m venv .venv
source .venv/bin/activate
python -m pip install -r requirements.lock.txt
python -m ipykernel install --user --name csi-genome --display-name "CSI Genome (Python 3.12)"
jupyter notebook CSI_Genome_Project07.ipynb
```

Select the **CSI Genome (Python 3.12)** kernel. The notebook has all functions and computation cells.
Choose **Restart Kernel and Run All Cells** to run it from the beginning.

- **Full analysis:** leave the default mode as `full` in the first code cell. This regenerates `results/`.
- **Quick demonstration:** change the first cell's mode to `test`. This regenerates `results_test/` and uses the same methods offline.

The full run downloads about 1 GB of compressed genomes and small FASTQ prefixes on the first run.
Allow about 6 GB of working disk space. Cached analysis takes a few minutes; downloading and NCBI waiting times vary.
Existing raw files and exact-query BLAST responses are reused. The notebook still recalculates the analysis.

The test uses **13 real assemblies and small samples of seven real read runs** (about 16 MB).
It should report **5/13 flagged**, with two adapter assemblies and three foreign candidates, and print **PASS** in Task 7.
Its rate is not an archive estimate. Some review cases are unavailable in the small dataset and remain unresolved.

## The eight tasks

| Task | What the notebook does |
|---|---|
| 1 | Defines foreign DNA, vectors, adapters, index hopping and sample swaps |
| 2 | Downloads the fixed genome/read inputs, provides optional metadata retrieval and checks identifier links |
| 3 | Scans adapters and compares that initial signature by platform and release era |
| 4 | Calculates short-read QC, trims reads and checks cross-sample assignments |
| 5 | Analyses long reads and searches for native/foreign split alignments |
| 6 | Aligns contigs, scores contamination versus HGT candidates and identifies selected sequences with BLAST |
| 7 | Calculates all category/taxon/year rates, reconstructs review evidence and checks the detector |
| 8 | Selects up to ten curator candidates with coordinates and evidence |

Task 3 starts with adapter signatures because foreign alignment is introduced in Task 6.
Task 7 then compares all categories. Test mode returns only the curator candidates actually present.

## Run without clicking cells

Snakemake executes the **same notebook** in order. It supplies the mode through `CSI_GENOME_MODE`.
The notebook contains the computation; Snakefile only starts it and tracks outputs.

```sh
snakemake --cores 2
```

This runs the offline test. For the complete study:

```sh
snakemake --cores 2 --config mode=full
```

Use `--forceall` to rerun a completed workflow. Executed workflow copies are stored locally in `.notebook_runs/`.

Docker alternative (definition supplied; image not validated):

```sh
docker build -t csi-genome .
docker run --rm -v "$PWD":/project csi-genome
```

## Data and review decisions

`config/assemblies.csv` records all assembly versions, BioSamples, BioProjects, platforms and release years.
`config/reads.csv` lists the seven read runs. NCBI and ENA metadata snapshots and reference manifests are also supplied.
The fixed inputs were retrieved on 29 September 2026; the original study and BLAST checks were performed on 30 September.
A notebook rerun measures its own runtime; it does not create a new BLAST response when the exact query is cached.

Set `REFRESH_METADATA = True` in full mode to download current annotations separately into `data/metadata_refresh/`.
The fixed study metadata are preserved. Large raw data and temporary files are ignored by Git.

**Review verdicts require judgment.** `config/review_annotations.csv` holds the recorded case-by-case decisions and reasons.
Task 7 selects a seeded sample again, extracts sequence context, joins only matching coordinates and computes error estimates.
Unreviewed cases remain unresolved. To review a case yourself, inspect `review_evidence.csv` and edit the matching annotation.
The notebook does not automatically perform an independent human validation.

## Limits

- Only exact adapter motifs and four foreign-reference sources are covered; altered adapters and other contaminants can be missed.
- Foreign matches and unusual GC can be real HGT or intended constructs. Vector matches require separate review.
- Short and long reads are small prefixes from selected runs, not random samples from the archive.
- Cross-sample mapping does not uniquely identify index hopping without index reads and lane information.
- `nt` was requested, but NCBI returned `core_nt`, as recorded in the saved responses.
- Reviewed false labels estimate a conditional false-discovery fraction, not a conventional false-positive rate on all clean genomes.

## Team contribution

Abylaikhan Torekhan and Kuanysh Bakizhan both worked on the code and checked the results.
Both partners are responsible for the final report and explaining the full analysis at the defence.
