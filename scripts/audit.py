"""A small archive audit: scan FASTA, align candidates, then count flags.
Coordinates in output tables are 1-based and inclusive.
"""
import argparse, csv, gzip, json, math, re, resource, time, zipfile
from pathlib import Path

ADAPTERS = {'TruSeq':'AGATCGGAAGAGCACACGTCT', 'TruSeq_R2':'AGATCGGAAGAGCGTCGTGTA', 'Nextera':'CTGTCTCTTATACACATCT',
            'BGI':'AAGTCGGAGGCCAAGCGGTC'}

def fasta(text):
    name, parts = None, []
    for line in text.splitlines():
        if line.startswith('>'):
            if name: yield name, ''.join(parts).upper()
            name, parts = line[1:].split()[0], []
        else: parts.append(line.strip())
    if name: yield name, ''.join(parts).upper()

def genomes(raw):
    for path in sorted(Path(raw).glob('genomes/*.zip')):
        with zipfile.ZipFile(path) as z:
            for name in z.namelist():
                if name.endswith('.fna'):
                    yield name.split('/')[-2], list(fasta(z.read(name).decode()))

def write_csv(path, rows, fields=None):
    path = Path(path); path.parent.mkdir(parents=True, exist_ok=True)
    with path.open('w') as f:
        w = csv.DictWriter(f, fieldnames=fields or list(rows[0]))
        w.writeheader(); w.writerows(rows)

def read_csv(path):
    with open(path) as f: return list(csv.DictReader(f))

def gc(s):
    n = sum(s.count(x) for x in 'ACGT')
    return (s.count('G') + s.count('C')) / n if n else 0

def rc(s): return s.translate(str.maketrans('ACGT','TGCA'))[::-1]

def scan(raw, out, metadata):
    start = time.perf_counter(); out = Path(out); out.mkdir(parents=True, exist_ok=True)
    meta = {r['accession']:r for r in read_csv(metadata)}
    motifs = {seq:name for name,s in ADAPTERS.items() for seq in (s,rc(s))}
    regex = re.compile('|'.join(motifs))
    rows, flags, candidates = [], [], []
    with (out/'candidates.fasta').open('w') as cf:
        for accession, contigs in genomes(raw):
            total = sum(len(s) for _,s in contigs)
            genome_gc = sum((s.count('G')+s.count('C')) for _,s in contigs) / max(1,total)
            adapter_count = outliers = 0
            for cid,s in contigs:
                delta = abs(gc(s)-genome_gc); hits = list(regex.finditer(s))
                for h in hits:
                    adapter_count += 1
                    flags.append(dict(accession=accession,contig=cid,category='Adapter',start=h.start()+1,end=h.end(),
                                      length=len(s),identity=1.0,coverage=len(h.group())/len(s),source=motifs[h.group()],score=4,
                                      status='technical_match',gc_delta=delta))
                # Screen every contig >=1 kb against a small panel; GC is evidence, not a gate.
                if len(s)>=1000:
                    qid=f'q{len(candidates):06d}'
                    candidates.append(dict(qid=qid,accession=accession,contig=cid,length=len(s),gc_delta=delta))
                    cf.write(f'>{qid}\n{s}\n')
                    outliers += delta>=0.12
            rows.append(dict(**meta[accession],bp=total,contigs=len(contigs),gc=genome_gc,
                             adapter_matches=adapter_count,gc_outlier_contigs=outliers))
    write_csv(out/'assemblies.csv',rows)
    write_csv(out/'adapter_hits.csv',flags,['accession','contig','category','start','end','length','identity','coverage','source','score','status','gc_delta'])
    write_csv(out/'candidates.csv',candidates)
    (out/'scan_runtime.json').write_text(json.dumps({'seconds':time.perf_counter()-start,'assemblies':len(rows),
       'bp':sum(r['bp'] for r in rows),'contigs':sum(r['contigs'] for r in rows), 'peak_rss_os_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},indent=2))
    print('Scanned',len(rows),'assemblies;',len(flags),'exact adapter matches',flush=True)

def align(out, panel, vectors):
    import mappy
    start=time.perf_counter(); out=Path(out)
    panel_index=mappy.Aligner(str(panel),preset='asm5',best_n=3)
    vector_index=mappy.Aligner(str(vectors),preset='asm5',best_n=3)
    if not panel_index or not vector_index: raise RuntimeError('Reference index failed')
    meta={r['qid']:r for r in read_csv(out/'candidates.csv')}
    assemblies={r['accession']:r for r in read_csv(out/'assemblies.csv')}
    rows=[]; sequences=[]
    for qid,seq,_ in mappy.fastx_read(str(out/'candidates.fasta')):
        r=meta[qid]; delta=float(r['gc_delta'])
        for category,index in [('Foreign',panel_index),('Vector',vector_index)]:
            hits=[h for h in index.map(seq) if h.mlen>=200 and h.mlen/max(1,h.blen)>=0.95 and h.mapq>=20]
            if not hits: continue
            # A circular chromosome may map in two pieces after rotation.
            by_target={}
            for hit in hits:by_target.setdefault(hit.ctg,[]).append(hit)
            group=max(by_target.values(),key=lambda hs:sum(h.mlen for h in hs))
            h=max(group,key=lambda h:h.mlen); intervals=sorted((h.q_st,h.q_en) for h in group)
            covered=0;stop=0
            for st,en in intervals:covered+=max(0,en-max(st,stop));stop=max(stop,en)
            left=min(st for st,en in intervals);right=max(en for st,en in intervals)
            ident=sum(h.mlen for h in group)/sum(h.blen for h in group);cov=covered/len(seq)
            flanks=left>=500 and len(seq)-right>=500
            score=2*(cov>=0.8)+(ident>=0.98)+(delta>=0.12)-3*flanks
            status='suspect_contamination' if score>=3 else 'HGT_or_construct_candidate' if flanks else 'uncertain'
            if category=='Vector': status='vector_match_needs_review'
            rows.append(dict(accession=r['accession'],contig=r['contig'],category=category,start=left+1,end=right,
                             length=len(seq),identity=ident,coverage=cov,source=h.ctg,score=int(score),status=status,gc_delta=delta))
            if category=='Foreign' and status=='suspect_contamination':
                mid=(h.q_st+h.q_en)//2; s=max(h.q_st,mid-500); e=min(h.q_en,s+1000)
                sequences.append((r['accession'],r['contig'],seq[s:e],s+1,e,len(seq)))
    write_csv(out/'alignment_hits.csv',rows,['accession','contig','category','start','end','length','identity','coverage','source','score','status','gc_delta'])
    # The 12 longest suspicious contigs, one per assembly; fixed deterministic order.
    sequences.sort(key=lambda r:(-r[5],r[0],r[1])); selected=[]; used=set()
    with (out/'nt_queries.fasta').open('w') as f:
        for acc,cid,s,st,en,n in sequences:
            if acc in used:continue
            used.add(acc); qid=f'Q{len(selected)+1:02d}'
            selected.append(dict(query=qid,accession=acc,contig=cid,start=st,end=en,length=n))
            f.write(f'>{qid}\n{s}\n')
            if len(selected)==12:break
    write_csv(out/'nt_queries.csv',selected)
    elapsed=time.perf_counter()-start
    (out/'align_runtime.json').write_text(json.dumps({'seconds':elapsed,'candidates':len(meta),
        'peak_rss_os_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},indent=2))
    print('Aligned',len(meta),'contigs;',len(rows),'panel/vector matches;',round(elapsed,1),'seconds',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('step',choices=['scan','align']);p.add_argument('--raw',default='data/raw')
    p.add_argument('--out',default='results');p.add_argument('--metadata',default='config/assemblies.csv')
    p.add_argument('--panel',default='data/raw/refs/panel.fasta');p.add_argument('--vectors',default='data/raw/refs/UniVec_Core.fasta')
    a=p.parse_args()
    if a.step=='scan':scan(a.raw,a.out,a.metadata)
    else:align(a.out,a.panel,a.vectors)
