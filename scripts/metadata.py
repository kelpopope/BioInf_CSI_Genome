"""Re-fetch the versioned records and map archive identifiers (optional live refresh)."""
import argparse,json,time,urllib.parse
from pathlib import Path
from fetch import download
from audit import read_csv,write_csv

def main():
    p=argparse.ArgumentParser();p.add_argument('--out',default='data/metadata_refresh');a=p.parse_args();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    accs=[r['accession'] for r in read_csv('config/assemblies.csv')]
    for i in range(0,len(accs),100):
        url='https://api.ncbi.nlm.nih.gov/datasets/v2/genome/accession/'+','.join(accs[i:i+100])+'/dataset_report?page_size=1000'
        download(url,out/f'assemblies_{i//100:03d}.json');time.sleep(.4)
    fields='run_accession,sample_accession,secondary_sample_accession,study_accession,instrument_platform,instrument_model,read_count,base_count,fastq_ftp,first_public'
    url='https://www.ebi.ac.uk/ena/portal/api/search?'+urllib.parse.urlencode(dict(result='read_run',query='tax_tree(31969)',fields=fields,format='tsv',limit=0))
    download(url,out/'ena_runs.tsv')
    # Saved config is the fixed 2026-09-29 snapshot. Live refresh is kept separate because metadata changes.
    print('Refreshed metadata in',out)

if __name__=='__main__':main()
