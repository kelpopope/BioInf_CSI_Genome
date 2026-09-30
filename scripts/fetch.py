"""Fetch the fixed accession list. Large files are cached outside Git."""
import argparse, csv, gzip, hashlib, io, json, shutil, time, urllib.parse, urllib.request
from pathlib import Path
from audit import fasta, read_csv

BASE='https://api.ncbi.nlm.nih.gov/datasets/v2/genome/download'

def download(url,path,payload=None,records=None):
    path=Path(path); path.parent.mkdir(parents=True,exist_ok=True)
    if path.exists() and path.stat().st_size:return
    request=urllib.request.Request(url,data=json.dumps(payload).encode() if payload else None,
        headers={'Content-Type':'application/json','User-Agent':'CSI-Genome-course-project'})
    temp=path.with_suffix(path.suffix+'.part')
    try:
        with urllib.request.urlopen(request,timeout=180) as response:
            if records:
                with gzip.GzipFile(fileobj=response) as source, gzip.open(temp,'wb') as out:
                    count=0
                    for _ in range(records):
                        lines=[source.readline() for _ in range(4)]
                        if not lines[0]:break
                        if not lines[0].startswith(b'@') or not lines[2].startswith(b'+') or len(lines[1].strip())!=len(lines[3].strip()):
                            raise ValueError('Incomplete or invalid FASTQ record')
                        out.writelines(lines); count+=1
                    if not count:raise ValueError('Empty FASTQ')
            else:
                with temp.open('wb') as out:shutil.copyfileobj(response,out)
        temp.replace(path)
    except Exception:
        temp.unlink(missing_ok=True); raise
    with open(path.parent/'downloads.jsonl','a') as log:
        log.write(json.dumps({'url':url,'file':path.name,'retrieved_utc':time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime()),
                             'sha256':hashlib.sha256(path.read_bytes()).hexdigest(),'first_records':records})+'\n')

def main():
    p=argparse.ArgumentParser();p.add_argument('--raw',default='data/raw');a=p.parse_args();raw=Path(a.raw)
    accessions=[r['accession'] for r in read_csv('config/assemblies.csv')]
    for i in range(0,len(accessions),150):
        download(BASE,raw/f'genomes/batch_{i//150:03d}.zip',{'accessions':accessions[i:i+150],'include_annotation_type':['GENOME_FASTA']})
    for r in read_csv('config/reads.csv'):
        for u in r['fastq_ftp'].split(';'):
            if r['library_layout']=='PAIRED' and not u.endswith(('_1.fastq.gz','_2.fastq.gz')):continue
            download('https://'+u,raw/'reads'/r['run_accession']/u.split('/')[-1],records=2000 if r['purpose']=='long_read' else 10000)
    # Use the exact small reference snapshot bundled with the real test dataset.
    (raw/'refs').mkdir(parents=True,exist_ok=True)
    for name in ['UniVec_Core.fasta','panel.fasta']:
        frozen=Path('data/test/refs')/name
        if frozen.exists() and not (raw/'refs'/name).exists():shutil.copyfile(frozen,raw/'refs'/name)
    download('https://ftp.ncbi.nlm.nih.gov/pub/UniVec/UniVec_Core',raw/'refs/UniVec_Core.fasta')
    if not (raw/'refs/panel.fasta').exists():
        import zipfile
        ids={'GCF_000005845.2':'Escherichia_coli','GCF_902374465.2':'Ralstonia_pickettii'}
        download(BASE,raw/'refs/panel.zip',{'accessions':list(ids),'include_annotation_type':['GENOME_FASTA']})
        with open(raw/'refs/panel.fasta','w') as out:
            with zipfile.ZipFile(raw/'refs/panel.zip') as z:
                for name in z.namelist():
                    if name.endswith('.fna'):
                        acc=name.split('/')[-2]
                        for cid,s in fasta(z.read(name).decode()):out.write(f'>{acc}|{ids[acc]}|{cid}\n{s}\n')
            for acc,label in [('NC_005089.1','Mus_musculus_mitochondrion'),('NC_000932.1','Arabidopsis_chloroplast')]:
                url='https://eutils.ncbi.nlm.nih.gov/entrez/eutils/efetch.fcgi?'+urllib.parse.urlencode(dict(db='nuccore',id=acc,rettype='fasta',retmode='text',email='abylai0209@gmail.com',tool='csi_genome'))
                path=raw/'refs'/f'{acc}.fasta';download(url,path)
                for cid,s in fasta(path.read_text()):out.write(f'>{acc}|{label}|{cid}\n{s}\n')
                time.sleep(.4)
    (raw/'.ready').touch()

if __name__=='__main__':main()
