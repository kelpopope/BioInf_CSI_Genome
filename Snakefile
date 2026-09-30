# Run exactly the same notebook, with a small offline dataset by default.
from pathlib import Path

MODE = config.get("mode", "test")
assert MODE in {"test", "full"}, "mode must be test or full"
OUT = "results_test" if MODE == "test" else "results"
RUN_NOTEBOOK = str(Path(f".notebook_runs/{MODE}.ipynb").resolve())

rule all:
    input:
        f"{OUT}/summary.json",
        f"{OUT}/nt_hits.csv",
        f"{OUT}/curator_top10.csv",
        f"{OUT}/review_summary.json",
        RUN_NOTEBOOK

rule analyse:
    input:
        notebook="CSI_Genome_Project07.ipynb",
        records="config/assemblies.csv",
        runs="config/reads.csv",
        review="config/review_annotations.csv",
        test_records="data/test/assemblies.csv",
        test_runs="data/test/reads.csv",
        references=["data/test/refs/panel.fasta", "data/test/refs/UniVec_Core.fasta"]
    output:
        summary=f"{OUT}/summary.json",
        blast=f"{OUT}/nt_hits.csv",
        curator=f"{OUT}/curator_top10.csv",
        review=f"{OUT}/review_summary.json",
        notebook=RUN_NOTEBOOK
    shell:
        "mkdir -p .notebook_runs && "
        "CSI_GENOME_MODE={MODE} python -m nbconvert --to notebook --execute "
        "{input.notebook:q} --output {output.notebook:q} "
        "--ExecutePreprocessor.timeout=1800"
