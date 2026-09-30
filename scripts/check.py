"""Check the real mini-dataset and measure an intentionally simple adapter baseline."""
import json,re
from pathlib import Path
from audit import ADAPTERS,genomes,rc

def main():
    result=json.loads(Path('results_test/summary.json').read_text())
    expected=json.loads(Path('data/test/expected.json').read_text())
    for key,value in expected.items():assert result[key]==value,(key,result[key],value)
    pattern=re.compile('|'.join(s for motif in ADAPTERS.values() for s in (motif,rc(motif))))
    sequence=next(cs[0][1] for acc,cs in genomes('data/test') if acc=='GCA_000008305.1')
    control=perfect=mutated=0;motif=ADAPTERS['Nextera']
    mutation=motif[:9]+('A' if motif[9]!='A' else 'C')+motif[10:]
    for i in range(100):
        s=sequence[1000*i:1000*(i+1)]
        control+=bool(pattern.search(s));perfect+=bool(pattern.search(s+motif));mutated+=bool(pattern.search(s+mutation))
    observed=dict(real_native_windows=100,native_flags=control,perfect_adapter_spikes_detected=perfect,
                  single_mismatch_spikes_detected=mutated,spikes_per_condition=100)
    assert control==0 and perfect==100
    Path('results/adapter_benchmark.json').write_text(json.dumps(observed,indent=2))
    print('Real test dataset: PASS. Adapter benchmark:',observed)

if __name__=='__main__':main()
