# Quick: snakemake --cores 2
# Full:  snakemake --cores 2 --config mode=full
MODE = config.get('mode', 'test')
RAW = 'data/test' if MODE == 'test' else 'data/raw'
OUT = 'results_test' if MODE == 'test' else 'results'
META = 'data/test/assemblies.csv' if MODE == 'test' else 'config/assemblies.csv'
RUNS = 'data/test/reads.csv' if MODE == 'test' else 'config/reads.csv'
OFFLINE = '--offline' if MODE == 'test' else ''
N_SHORT = 100 if MODE == 'test' else 10000
N_LONG = 20 if MODE == 'test' else 2000

rule all:
    input: f'{OUT}/summary.json', f'{OUT}/nt_hits.csv', f'{OUT}/curator_top10.csv'

rule fetch:
    output: 'data/raw/.ready'
    shell: 'python scripts/fetch.py'

rule scan:
    input: ready=f'{RAW}/.ready', metadata=META, code='scripts/audit.py'
    output: f'{OUT}/assemblies.csv', f'{OUT}/adapter_hits.csv', f'{OUT}/candidates.csv', f'{OUT}/candidates.fasta'
    shell: 'python scripts/audit.py scan --raw {RAW} --out {OUT} --metadata {META}'

rule align:
    input: f'{OUT}/candidates.csv', f'{OUT}/candidates.fasta', f'{OUT}/assemblies.csv', 'scripts/audit.py'
    output: f'{OUT}/alignment_hits.csv', f'{OUT}/nt_queries.fasta', f'{OUT}/nt_queries.csv'
    shell: 'python scripts/audit.py align --out {OUT} --panel {RAW}/refs/panel.fasta --vectors {RAW}/refs/UniVec_Core.fasta'

rule reads:
    input: ready=f'{RAW}/.ready', runs=RUNS, code='scripts/reads.py', common='scripts/audit.py'
    output: f'{OUT}/read_qc.csv', f'{OUT}/read_gc.csv', f'{OUT}/read_quality.csv', f'{OUT}/index_hopping.csv', f'{OUT}/chimeric_reads.csv'
    shell: 'python scripts/reads.py --raw {RAW} --out {OUT} --runs {RUNS} --panel {RAW}/refs/panel.fasta --short {N_SHORT} --long {N_LONG}'

rule blast:
    input: f'{OUT}/nt_queries.fasta', f'{OUT}/nt_queries.csv', 'scripts/blast_nt.py'
    output: f'{OUT}/nt_hits.csv'
    shell: 'python scripts/blast_nt.py --out {OUT} {OFFLINE}'

rule summary:
    input: f'{OUT}/assemblies.csv', f'{OUT}/adapter_hits.csv', f'{OUT}/alignment_hits.csv',
           f'{OUT}/read_qc.csv', f'{OUT}/nt_hits.csv', 'scripts/summarize.py'
    output: f'{OUT}/summary.json', f'{OUT}/curator_top10.csv'
    shell: 'python scripts/summarize.py --out {OUT}'
