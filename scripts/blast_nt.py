"""Scripted NCBI nt BLAST with an exact-query cache and bounded waiting."""
import argparse, hashlib, json, re, time, urllib.parse, urllib.request
from pathlib import Path
from audit import read_csv, write_csv
URL='https://blast.ncbi.nlm.nih.gov/Blast.cgi'

def request(params):
    data=urllib.parse.urlencode(params).encode()
    with urllib.request.urlopen(urllib.request.Request(URL,data=data),timeout=180) as r:return r.read().decode()

def main(a):
    out=Path(a.out);query=(out/'nt_queries.fasta').read_text();digest=hashlib.sha256(query.encode()).hexdigest()
    cache=Path(a.cache);cache.mkdir(parents=True,exist_ok=True);result=cache/f'{digest}.json';stamp=cache/f'{digest}.request.json'
    if not result.exists():
        if a.offline:raise FileNotFoundError('No BLAST response for these exact queries')
        if stamp.exists():rid=json.loads(stamp.read_text())['rid']
        else:
            params=dict(CMD='Put',PROGRAM='blastn',DATABASE='nt',QUERY=query,MEGABLAST='on',HITLIST_SIZE=20,EXPECT='1e-10',
                        tool='csi_genome_student_project',email='abylai0209@gmail.com')
            text=request(params);match=re.search(r'RID = (\S+)',text)
            if not match:raise RuntimeError('NCBI did not return a BLAST request ID')
            rid=match.group(1)
            stamp.write_text(json.dumps(dict(rid=rid,query_sha256=digest,database_requested='nt',program='megablast',
                 submitted_utc=time.strftime('%Y-%m-%dT%H:%M:%SZ',time.gmtime())),indent=2))
        print('BLAST request:',rid,flush=True)
        for _ in range(30):
            time.sleep(60)
            status=request(dict(CMD='Get',FORMAT_OBJECT='SearchInfo',RID=rid))
            if 'Status=READY' in status:
                data=json.loads(request(dict(CMD='Get',FORMAT_TYPE='JSON2_S',RID=rid)))
                result.write_text(json.dumps(data));break
            if 'Status=FAILED' in status or 'Status=UNKNOWN' in status:raise RuntimeError('BLAST failed/expired: '+rid)
        else:raise TimeoutError('BLAST still queued. Re-run to resume the same request.')
    data=json.loads(result.read_text());queries={r['query']:r for r in read_csv(out/'nt_queries.csv')};rows=[]
    for item in data['BlastOutput2']:
        report=item['report'];s=report['results']['search'];qid=s['query_title'].split()[0];own=queries[qid]['contig'].split('.')[0]
        best=[]
        for hit in s.get('hits',[]):
            desc=next((d for d in hit['description'] if d.get('accession','').split('.')[0]!=own),None)
            if desc is None:continue
            h=max(hit['hsps'],key=lambda h:h['bit_score'])
            best.append(dict(**queries[qid],subject=desc.get('accession',''),organism=desc.get('sciname',''),
                    title=desc.get('title',''),identity=h['identity']/h['align_len'],query_coverage=(h['query_to']-h['query_from']+1)/s['query_len'],
                    evalue=h['evalue'],bitscore=h['bit_score'],database=report['search_target'].get('db',''),query_sha256=digest))
        if best:rows.append(max(best,key=lambda r:r['bitscore']))
        else:rows.append(dict(**queries[qid],subject='',organism='No non-self hit',title='',identity=0,query_coverage=0,evalue='',bitscore=0,
                              database=report['search_target'].get('db',''),query_sha256=digest))
    write_csv(out/'nt_hits.csv',rows)
    print('Saved',len(rows),'BLAST results',flush=True)

if __name__=='__main__':
    p=argparse.ArgumentParser();p.add_argument('--out',default='results');p.add_argument('--cache',default='data/blast_cache');p.add_argument('--offline',action='store_true');main(p.parse_args())
