"""QC for small real FASTQ prefixes; mapping suggests artefacts, not their cause."""
import resource, argparse, collections, gzip, itertools, json, math, time
from pathlib import Path
import mappy, numpy as np
from audit import ADAPTERS, fasta, genomes, gc, rc, read_csv, write_csv

def fastq(path,n):
    with gzip.open(path,'rt') as f:
        for _ in range(n):
            name=f.readline().strip()
            if not name:return
            seq=f.readline().strip();plus=f.readline();qual=f.readline().strip()
            if not name.startswith('@') or not plus.startswith('+') or len(seq)!=len(qual):raise ValueError('Invalid FASTQ')
            yield name[1:].split()[0],seq,qual

def entropy(s):
    counts=collections.Counter(s[i:i+3] for i in range(len(s)-2)); n=sum(counts.values())
    return -sum(c/n*math.log2(c/n) for c in counts.values()) if n else 0

def main(a):
    start=time.perf_counter();out=Path(a.out);out.mkdir(parents=True,exist_ok=True)
    runs=read_csv(a.runs); wanted={r['accession'] for r in runs};own={}
    for acc,contigs in genomes(a.raw):
        if acc in wanted:own[acc]=contigs
    panel=Path(a.panel).read_text();summary=[];hist=[];cycle=[];splits=[];pool=[]
    # Cross-sample index combines two candidate batch mates; release date is only a proxy for pooling.
    pairs=[r for r in runs if r['purpose']=='index_hopping']
    poolfa=out/'pooled.fasta'
    with poolfa.open('w') as f:
        for r in pairs:
            for cid,s in own[r['accession']]:f.write(f'>{r["run_accession"]}|{cid}\n{s}\n')
    poolindex=mappy.Aligner(str(poolfa),preset='sr') if pairs else None
    for r in runs:
        run=r['run_accession'];long=r['purpose']=='long_read';n=a.long if long else a.short
        preset='map-pb' if r['platform']=='PacBio' and 'RS II' in r['instrument_model'] else 'map-hifi' if r['platform']=='PacBio' else 'map-ont' if long else 'sr'
        ref=out/f'{run}.reference.fasta'
        with ref.open('w') as f:
            for cid,s in own[r['accession']]:f.write(f'>native|{cid}\n{s}\n')
            f.write(panel)
        index=mappy.Aligner(str(ref),preset=preset,best_n=5)
        paths=sorted((Path(a.raw)/'reads'/run).glob('*.fastq.gz'))
        paths=[p for p in paths if p.name.endswith(('_1.fastq.gz','_2.fastq.gz'))] if not long else paths[:1]
        lengths=[];qualities=[];gcs=[];identities=[];indels=[];seen=set();duplicates=adapter=low=polyg=mapped=chim=0
        trimmed=discarded=0; retained_q=[]
        qs=collections.defaultdict(list);cross=paircount=0
        for file in paths:
            for name,s,q in fastq(file,n):
                lengths.append(len(s));gcs.append(gc(s));qualities.append(np.mean([ord(x)-33 for x in q]));low+=entropy(s)<3
                duplicates+=s in seen;seen.add(s);polyg+=s.endswith('G'*10)
                adapter+=any(x in s for x in ADAPTERS.values())
                for pos,char in enumerate(q):
                    if not long:qs[pos].append(ord(char)-33)
                if not long:
                    cut=min([s.find(x) for x in ADAPTERS.values() if x in s] or [len(s)])
                    clean=s[:cut]; quality=q[:cut]
                    while quality and ord(quality[-1])-33<20:clean=clean[:-1];quality=quality[:-1]
                    trimmed+=len(clean)<len(s)
                    if len(clean)<50:discarded+=1;continue
                    retained_q.extend(ord(x)-33 for x in quality)
                    s=clean
                hits=list(index.map(s,cs=True));good=[h for h in hits if h.mapq>=20 and h.blen>=100]
                if good:
                    h=max(good,key=lambda h:h.mlen);mapped+=1;identities.append(h.mlen/h.blen)
                    insdel=sum(k for k,op in (h.cigar or []) if op in [1,2]);errors=h.blen-h.mlen
                    if errors:indels.append(insdel/errors)
                if long:
                    found=False
                    for h1,h2 in itertools.combinations(good,2):
                        overlap=max(0,min(h1.q_en,h2.q_en)-max(h1.q_st,h2.q_st))
                        different=h1.ctg.startswith('native|')!=h2.ctg.startswith('native|')
                        if different and min(h1.q_en-h1.q_st,h2.q_en-h2.q_st)>=300 and overlap<=50:
                            chim+=1;found=True
                            splits.append(dict(run=run,read=name,left=h1.ctg,right=h2.ctg,left_start=h1.q_st+1,left_end=h1.q_en,right_start=h2.q_st+1,right_end=h2.q_en))
                            break
        if r['purpose']=='index_hopping' and len(paths)==2:
            for left,right in zip(fastq(paths[0],n),fastq(paths[1],n),strict=True):
                paircount+=1;assigned=[]
                for _,s,_ in [left,right]:
                    hs=[h for h in poolindex.map(s) if h.mapq>=30 and h.mlen/max(1,h.blen)>=.98 and (h.q_en-h.q_st)>=.9*len(s)]
                    assigned.append(hs[0].ctg.split('|')[0] if hs else '')
                cross+=bool(assigned[0] and assigned[0]==assigned[1] and assigned[0]!=run)
            pool.append(dict(run=run,pairs=paircount,cross_sample_pairs=cross,fraction=cross/max(1,paircount)))
        size=len(lengths)
        summary.append(dict(run=run,assembly=r['accession'],platform=r['platform'],model=r['instrument_model'],preset=preset,reads=size,
          median_length=float(np.median(lengths)),mean_phred=float(np.mean(qualities)),adapter_fraction=adapter/size,low_complexity_fraction=low/size,
          duplicate_fraction=duplicates/size,polyg_fraction=polyg/size,mapped_fraction=mapped/size,
          trimmed_fraction=trimmed/size,discarded_fraction=discarded/size,mean_retained_phred=float(np.mean(retained_q)) if retained_q else None,
          median_identity=float(np.median(identities)) if identities else None,
          median_indel_error_share=float(np.median(indels)) if indels else None,foreign_native_splits=chim,chimera_fraction=chim/size))
        for pos,values in qs.items():cycle.append(dict(run=run,cycle=pos+1,mean_phred=float(np.mean(values))))
        counts,edges=np.histogram(gcs,bins=np.linspace(0,1,21))
        for i,c in enumerate(counts):hist.append(dict(run=run,gc_mid=(edges[i]+edges[i+1])/2,count=int(c)))
        print(run,size,'reads; split candidates:',chim,flush=True)
    write_csv(out/'read_qc.csv',summary);write_csv(out/'read_gc.csv',hist);write_csv(out/'read_quality.csv',cycle)
    write_csv(out/'index_hopping.csv',pool,['run','pairs','cross_sample_pairs','fraction'])
    write_csv(out/'chimeric_reads.csv',splits,['run','read','left','right','left_start','left_end','right_start','right_end'])
    (out/'reads_runtime.json').write_text(json.dumps({'seconds':time.perf_counter()-start,'peak_rss_os_units':resource.getrusage(resource.RUSAGE_SELF).ru_maxrss},indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--raw',default='data/raw');p.add_argument('--out',default='results')
    p.add_argument('--runs',default='config/reads.csv');p.add_argument('--panel',default='data/raw/refs/panel.fasta')
    p.add_argument('--short',type=int,default=10000);p.add_argument('--long',type=int,default=2000)
    main(p.parse_args())
