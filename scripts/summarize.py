"""Create the tables and figures used in the notebook and report."""
import argparse, json, math, random
from pathlib import Path
import pandas as pd
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
from audit import read_csv, write_csv

COLORS=['#20639b','#ed553b','#3caea3','#7c6a9c']
plt.rcParams.update({'font.size':10,'axes.spines.top':False,'axes.spines.right':False,'figure.dpi':160})

def interval(k,n):
    if not n:return [None,None]
    z=1.96;p=k/n;den=1+z*z/n;mid=(p+z*z/(2*n))/den;h=z*math.sqrt(p*(1-p)/n+z*z/(4*n*n))/den
    return [max(0,mid-h),min(1,mid+h)]

def main(a):
    out=Path(a.out);fig=out/'figures';fig.mkdir(exist_ok=True)
    A=pd.read_csv(out/'assemblies.csv');D=pd.read_csv(out/'adapter_hits.csv');V=pd.read_csv(out/'alignment_hits.csv')
    F=V[(V.category=='Foreign') & (V.status=='suspect_contamination')]
    A['adapter']=A.accession.isin(D.accession)
    A['foreign']=A.accession.isin(F.accession)
    A['vector']=A.accession.isin(V[V.category=='Vector'].accession)
    A['hgt']=A.accession.isin(V[V.status=='HGT_or_construct_candidate'].accession)
    A['flagged']=A.adapter|A.foreign
    A['taxon']=A.organism.str.replace(r'^(uncultured |Candidatus )+','',regex=True).str.split().str[0].str.strip("'")
    A['era']=pd.cut(A.year,[1990,2015,2019,2022,2026],labels=['<=2015','2016-2019','2020-2022','2023-2026'])
    A.to_csv(out/'assembly_flags.csv',index=False)
    rate_rows=[]
    for by in ['year','era','platform','taxon']:
        for label,g in A.groupby(by,observed=True):
            for cat in ['adapter','foreign','vector','hgt','flagged']:
                k=int(g[cat].sum());n=len(g);lo,hi=interval(k,n)
                rate_rows.append(dict(grouping=by,group=label,category=cat,n=n,flagged=k,rate=k/n,ci_low=lo,ci_high=hi))
    write_csv(out/'rates.csv',rate_rows)
    # Explicit identifier reconciliation against both ENA BioSample columns.
    ena=pd.read_csv('config/ena_runs.tsv',sep='\t',dtype=str).fillna('')
    links=[]
    for row in A.itertuples():
        if pd.isna(row.biosample):continue
        found=ena[(ena.sample_accession==row.biosample)|(ena.secondary_sample_accession==row.biosample)]
        for r in found.itertuples():links.append(dict(assembly=row.accession,biosample=row.biosample,run=r.run_accession,
            ena_sample=r.sample_accession,ena_secondary_sample=r.secondary_sample_accession,platform=r.instrument_platform))
    write_csv(out/'archive_links.csv',links,['assembly','biosample','run','ena_sample','ena_secondary_sample','platform'])
    n=len(A);k=int(A.flagged.sum())
    summary=dict(assemblies=n,contigs=int(A.contigs.sum()),bases=int(A.bp.sum()),adapter_assemblies=int(A.adapter.sum()),
        foreign_candidates=int(A.foreign.sum()),vector_candidates=int(A.vector.sum()),hgt_or_construct_candidates=int(A.hgt.sum()),
        flagged_assemblies=k,flag_rate=k/n,flag_rate_ci95=interval(k,n),ani_failed=int((A.ani_status=='Failed').sum()),
        direct_ena_links=len({r['assembly'] for r in links}),missing_platform=int((A.platform=='Unknown').sum()),
        definition='Flagged = exact known adapter motif OR high-scoring foreign panel match. This is not proven contamination prevalence.')
    for name in ['scan','align','reads']:
        p=out/f'{name}_runtime.json'
        if p.exists():summary[name+'_runtime']=json.loads(p.read_text())
    (out/'summary.json').write_text(json.dumps(summary,indent=2))
    # Curator list: foreign calls plus adapter-rich assemblies, each with a concrete coordinate.
    top=[]
    for row in F.sort_values('length',ascending=False).drop_duplicates('accession').itertuples():
        top.append(dict(accession=row.accession,category='Foreign',contig=row.contig,start=row.start,end=row.end,
            evidence=f'{100*row.identity:.2f}% identity; {100*row.coverage:.1f}% contig covered; {row.source}',
            action='Inspect BLAST, raw-read support and assembly context; remove only if confirmed.'))
    for acc,count in D.groupby('accession').size().sort_values(ascending=False).items():
        if len(top)>=10:break
        if acc in {r['accession'] for r in top}:continue
        row=D[D.accession==acc].iloc[0]
        top.append(dict(accession=acc,category='Adapter',contig=row.contig,start=row.start,end=row.end,
           evidence=f'{count} exact adapter matches; example: {row.source}',action='Check listed coordinates; trim terminal adapter or inspect internal join.'))
    write_csv(out/'curator_top10.csv',top)
    # Fixed random, assembly-level review sample; decisions are supplied separately after inspection.
    if not (out/'review_sample.csv').exists():
        selected=[]
        for name,frame,limit in [('Adapter',D,20),('Vector',V[V.category=='Vector'],10),('Foreign',F,10)]:
            unique=frame.sort_values(['accession','contig','start']).drop_duplicates('accession')
            if len(unique):selected.extend(unique.sample(n=min(limit,len(unique)),random_state=42).to_dict('records'))
        write_csv(out/'review_sample.csv',selected)
    # Figures are generated from the same tables as the text.
    fig1,axs=plt.subplots(1,2,figsize=(9,3.3))
    labels=['Adapter','Foreign','Vector*','HGT/construct*'];counts=[summary[x] for x in ['adapter_assemblies','foreign_candidates','vector_candidates','hgt_or_construct_candidates']]
    axs[0].bar(labels,counts,color=COLORS);axs[0].set_ylabel('Assemblies');axs[0].set_title('Screening categories (overlap allowed)')
    rates=pd.DataFrame(rate_rows);p=rates[(rates.grouping=='platform')&(rates.category=='flagged')]
    axs[1].bar(p.group,100*p.rate,color=COLORS[0]);axs[1].set_ylabel('Assemblies flagged (%)');axs[1].set_title('By declared sequencing platform')
    fig1.tight_layout();fig1.savefig(fig/'archive.png');plt.close(fig1)
    p=rates[(rates.grouping=='year')&(rates.category=='flagged')]
    f,ax=plt.subplots(figsize=(9,2.8));ax.bar(p.group.astype(int),100*p.rate,color=COLORS[0]);ax.set(xlabel='Assembly release year',ylabel='Flag rate (%)');f.tight_layout();f.savefig(fig/'years.png');plt.close(f)
    if (out/'read_qc.csv').exists():
        Q=pd.read_csv(out/'read_qc.csv');G=pd.read_csv(out/'read_gc.csv');B=pd.read_csv(out/'read_quality.csv')
        f,axs=plt.subplots(1,2,figsize=(9,3.2))
        for run,g in B.groupby('run'):axs[0].plot(g.cycle,g.mean_phred,label=run)
        axs[0].set(xlabel='Read cycle (R1 and R2 combined)',ylabel='Mean Phred score');axs[0].legend(fontsize=7)
        for run,g in G.groupby('run'):
            if run in ['SRR9681681','SRR33330896']:axs[1].plot(100*g.gc_mid,g['count']/g['count'].sum(),label=run)
        axs[1].set(xlabel='GC per read (%)',ylabel='Fraction of reads');axs[1].legend(fontsize=7)
        f.tight_layout();f.savefig(fig/'short_reads.png');plt.close(f)
        f,axs=plt.subplots(1,2,figsize=(9,3))
        x=range(len(Q));axs[0].bar(x,Q.median_length,color=COLORS[0]);axs[0].set_yscale('log');axs[0].set_ylabel('Median length (bp, log scale)')
        axs[1].bar(x,100*(1-Q.median_identity),color=COLORS[1]);axs[1].set_ylabel('Median alignment difference (%)')
        for ax in axs:ax.set_xticks(list(x),Q.run,rotation=35,ha='right',fontsize=7)
        f.tight_layout();f.savefig(fig/'platforms.png');plt.close(f)
    print(json.dumps(summary,indent=2))

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',default='results');main(p.parse_args())
