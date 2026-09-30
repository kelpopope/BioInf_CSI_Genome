# CSI Genome: The Contamination Detective

**Project 07 · Introduction to Bioinformatics · Astana IT University**  
**Abylaikhan Torekhan & Kuanysh Bakizhan · BDA-2405**

We screened **3,006 real GenBank Mollicutes assemblies** using a small, readable pipeline.
**246 assemblies (8.18%) were flagged:** 243 contained known adapter motifs and three had strong
foreign-DNA matches. NCBI BLAST supported E. coli in one case and mouse mitochondrial DNA in two.
These are screening flags, not a measurement of all contamination in GenBank.

## Read the work

- [Notebook](CSI_Genome_Project07.ipynb): question, methods, results and reproduction command.
- [Notebook as HTML](notebook.html).
- [12-page report](report/CSI_Genome_Report.pdf), including Appendix A on AI usage.
- [Editable report in Word](report/CSI_Genome_Report.docx): source for the final PDF.
  The earlier report builder generates a separate draft and does not replace the edited report.
- [Ten curator candidates](results/curator_top10.csv), with coordinates and suggested checks.
- [Case review](results/review_decisions.csv), including unresolved cases and intentional constructs.

## Run from a clean machine

Use Python 3.12 and a C compiler for mappy (Linux build-essential / macOS Command Line Tools).

1. Clone this repository and open its folder.
2. `python3 -m venv .venv`
3. Activate it: `source .venv/bin/activate` (Windows: `.venv\Scripts\activate`).
4. `python -m pip install -r requirements.lock.txt`
5. **Quick, offline reproduction:** `snakemake --cores 2`
6. **Check the expected result:** `python scripts/check.py`
7. **Full headline result from raw inputs:** `snakemake --cores 2 --config mode=full`

The bundled test has 13 real assemblies and short prefixes of seven real sequencing runs
(about 16 MB). Expected result: **5/13 flagged**, two adapter assemblies and three foreign candidates.
It uses the same code and cached BLAST answers as the full run; its flag rate is not an archive estimate.
The BLAST cache is matched to the exact query hash. The complete full run was also executed through Snakemake.

The full run downloads about 1 GB of compressed genomes, references and small FASTQ prefixes.
Cached analysis takes a few minutes on the tested laptop; downloads and the NCBI queue vary.
Budget about **6 GB of free disk space**. The raw data and large intermediate FASTA are ignored by Git.
`results/summary.json` and the runtime JSON files contain the measured results and timings.

Docker alternative:

```sh
docker build -t csi-genome .
docker run --rm -v "$PWD":/project csi-genome
```

## Understand the code

| File | Job |
|---|---|
| `scripts/fetch.py` | Download the fixed accession list and limited FASTQ prefixes |
| `scripts/metadata.py` | Optional live metadata refresh, saved separately from the study snapshot |
| `scripts/audit.py` | Scan adapters, map to a small foreign/vector panel, score candidates |
| `scripts/reads.py` | Short-read QC, long-read splits and cross-sample read check |
| `scripts/blast_nt.py` | Scripted NCBI search with exact-query caching and self-hit exclusion |
| `scripts/summarize.py` | Assemble rates, archive links, figures and curator candidates |
| `scripts/check.py` | Verify the mini-dataset and run the adapter-spike baseline |
| `Snakefile` | Connect the steps in the correct order |

All assembly accessions, dates, declared platforms, BioSamples and BioProjects are in
`config/assemblies.csv`. Read-run IDs and URLs are in `config/reads.csv`. NCBI/ENA metadata snapshots
are supplied in `config/`; original public archive downloads date to **2026-09-29**.
This new analysis and BLAST verification were performed on **2026-09-30**.
`scripts/metadata.py` demonstrates repeatable retrieval without silently replacing the fixed snapshot.

`nt` was requested from NCBI; the returned database is **core_nt**, as recorded in the raw responses.
Some large eukaryotic chromosomes are excluded from core_nt, which limits host-DNA detection.
Source links and biological limitations are in the report. No dataset was generated to stand in for real data.

## Main limitations

- The foreign panel has four references. Missing organisms and damaged/partial adapters can be missed.
- A vector match can be natural bacterial sequence or an intentional construct; it is not included in the main count.
- HGT, sample swaps and index hopping cannot be proved from the available evidence alone.
- Reads are the first 10,000 pairs or 2,000 long reads per run, not random archive-wide samples.
- Review estimates are conditional false-discovery fractions, with unresolved cases reported separately.
  They are evidence judgments, not independent experimental truth.

## Team responsibilities and commits

Abylaikhan Torekhan and Kuanysh Bakizhan both worked on the code and checked the results.
Both partners are responsible for the final report and for explaining the full analysis at the defence.

Joint commits use GitHub's `Co-authored-by` trailer with the verified account address
`258616596+karp3n3@users.noreply.github.com`. The repository preserves the actual commit dates;
co-author metadata does not itself establish who completed a technical task.
