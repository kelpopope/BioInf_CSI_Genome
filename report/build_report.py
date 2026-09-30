"""Build an earlier draft from result tables; the final report is edited in CSI_Genome_Report.docx."""
import json,sys
from pathlib import Path
import pandas as pd
from reportlab.platypus import SimpleDocTemplate,Paragraph,Spacer,PageBreak,Table,TableStyle,Image
from reportlab.lib.styles import getSampleStyleSheet,ParagraphStyle
from reportlab.lib import colors
from reportlab.lib.enums import TA_CENTER
from reportlab.lib.pagesizes import A4
from xml.sax.saxutils import escape
ROOT=Path(__file__).resolve().parents[1];R=ROOT/'results';H=json.loads((R/'summary.json').read_text())
A=pd.read_csv(R/'assembly_flags.csv');Q=pd.read_csv(R/'read_qc.csv');rates=pd.read_csv(R/'rates.csv');N=pd.read_csv(R/'nt_hits.csv')
style=getSampleStyleSheet();style.add(ParagraphStyle(name='BodyText2',fontName='Helvetica',fontSize=10.2,leading=14.4,spaceAfter=8,textColor=colors.HexColor('#263445')))
style.add(ParagraphStyle(name='Small2',fontName='Helvetica',fontSize=8.1,leading=10.6,spaceAfter=6))
style['Title'].fontSize=27;style['Title'].leading=31;style['Title'].textColor=colors.HexColor('#153c5b')
style['Heading1'].fontSize=19;style['Heading1'].leading=23;style['Heading1'].spaceAfter=14
style['Heading2'].fontSize=12;style['Heading2'].spaceBefore=10
story=[]
def p(text,small=False):story.append(Paragraph(text,style['Small2' if small else 'BodyText2']))
def title(t):story.append(Paragraph(t,style['Heading1']))
def sub(t):story.append(Paragraph(t,style['Heading2']))
def table(headers,rows,widths=None):
    data=[[Paragraph(escape(str(v)),style['Small2']) for v in row] for row in [headers]+rows]
    t=Table(data,colWidths=widths,repeatRows=1,hAlign='LEFT');t.setStyle(TableStyle([('BACKGROUND',(0,0),(-1,0),colors.HexColor('#e6eef4')),('VALIGN',(0,0),(-1,-1),'TOP'),('LINEBELOW',(0,0),(-1,0),.6,colors.HexColor('#7796ad')),('LINEBELOW',(0,1),(-1,-1),.2,colors.HexColor('#d6dfe5')),('LEFTPADDING',(0,0),(-1,-1),6),('RIGHTPADDING',(0,0),(-1,-1),6),('TOPPADDING',(0,0),(-1,-1),5),('BOTTOMPADDING',(0,0),(-1,-1),4)]));story.append(t);story.append(Spacer(1,9))
def figure(file,caption,height=180):story.append(Image(str(R/'figures'/file),width=510,height=height,kind='proportional'));p(caption,True)
def page():story.append(PageBreak())
def pct(v):return f'{100*v:.2f}%'

story.append(Paragraph('CSI Genome',style['Title']));p('Project 07 - The Contamination Detective')
p('<b>Abylaikhan Torekhan &amp; Kuanysh Bakizhan</b><br/>BDA-2405 - Astana IT University<br/>Introduction to Bioinformatics - 30 September 2026')
story.append(Spacer(1,14));title('A small, reproducible audit of 3,006 genomes')
p(f'<b>Main finding:</b> {H["flagged_assemblies"]} of {H["assemblies"]:,} GenBank Mollicutes assemblies ({pct(H["flag_rate"])}) were flagged by a deliberately narrow screen. This is a <b>screening rate</b>, not an estimate of every form of contamination in the archive.')
table(['Evidence','Assemblies','Meaning'],[['Known adapter motif',H['adapter_assemblies'],'An exact 19-21 bp technical sequence'],['Foreign DNA candidate',H['foreign_candidates'],'Near-complete match to a small foreign reference panel'],['Vector match (separate)',H['vector_candidates'],'Needs biological context; not automatically contamination'],['HGT / construct context',H['hgt_or_construct_candidates'],'Internal foreign block with long flanks']], [160,65,285])
p('BLAST identified E. coli DNA in one assembly and mouse mitochondrial DNA in two others. The report supplies ten curator candidates with versioned accessions, coordinates, evidence and suggested checks. It does not recommend deleting sequence solely on the basis of a flag.')
p('The project uses a few readable Python scripts and a Snakemake workflow. The real mini-dataset runs offline. The full archive analysis uses cached downloads and small prefixes of seven sequencing runs. No neural network, distributed cluster or full local nt database is required.')
sub('What this audit adds')
p('It separates evidence of a foreign sequence from evidence that the sequence is an error. In a case review, several strong vector matches belonged to deliberately engineered strains. Counting every vector match as contamination would therefore be misleading.')
page()

title('1. Biological question and study design')
p('How often do assemblies in a defined public archive slice contain recognizable technical sequence or strong matches to unrelated organisms? Which records deserve curator attention, and which apparent problems are better explained by engineering or horizontal gene transfer?')
p('The fixed slice contains every versioned GenBank assembly returned for Mollicutes (NCBI taxid 31969) in the saved 29 September 2026 snapshot. The 3,006 records contain 317,568 contigs and 3,674,465,179 bases. Small, often AT-rich mollicute genomes make this a manageable teaching dataset. Many organisms are associated with animal or plant hosts, so host DNA is biologically plausible.')
table(['Category','Definition and laboratory origin'],[['Foreign organism','DNA from a host, another microbe, reagents or a mixed metagenomic bin.'],['Vector','Cloning plasmid or construct sequence. It can be accidental carry-over or an intentional engineered insertion.'],['Adapter','Library-preparation sequence retained after read-through, incomplete trimming or an assembly join.'],['Index hopping','A sequencing read receives an index belonging to another library during multiplex sequencing.'],['Sample swap / mislabel','The specimen, species label or archive link refers to the wrong organism.']], [110,400])
p('Horizontal gene transfer is a real biological movement of DNA between organisms. An internal foreign block with substantial flanking sequence is therefore treated as an HGT-or-construct candidate, not automatically removed. Flanks alone do not prove HGT: a chimeric assembly can produce the same pattern.')
sub('Scope decided before interpretation')
p('The assembly screen uses known adapter strings, UniVec_Core and four foreign references: E. coli, Ralstonia pickettii, mouse mitochondria and an Arabidopsis chloroplast. These represent laboratory, reagent, animal-host and plant-host sources. This small panel is fast and explainable but has limited sensitivity. A negative result means that this screen found no evidence, not that the assembly is clean.')
p('Real archived sequences are the main data. Additional synthetic adapter spikes are used only to measure the behaviour of the exact-match detector. They do not replace the archive audit.')
page()

title('2. Data sources and identifier mapping')
table(['Source','Data and identifiers'],[['GenBank / NCBI Datasets','3,006 versioned GCA accessions, FASTA packages, BioSample and BioProject metadata. Full list: config/assemblies.csv.'],['ENA','Read-run metadata and FASTQ. Links use both sample_accession and secondary_sample_accession.'],['RefSeq / Nucleotide','GCF_000005845.2, GCF_902374465.2, NC_005089.1 and NC_000932.1.'],['UniVec_Core','NCBI reference of vector, adapter and related sequences; saved reference is checksummed.'],['NCBI BLAST','Three 1 kb queries. Request asked for nt; the returned database was core_nt. Exact queries, raw JSON responses and request IDs are saved.']], [135,375])
p(f'Direct BioSample matching linked {H["direct_ena_links"]} assemblies to ENA runs. This deliberately simple join misses MAGs whose reads are registered under a parent metagenome sample. Missing links are retained as missing, not interpreted as absent sequencing. {H["missing_platform"]} assemblies had unrecognized or missing platform text under the simple normalization rule. Seven records have an NCBI ANI taxonomy status of Failed; they are taxonomy-review candidates, not seven demonstrated tube swaps.')
table(['Run','Platform / model','Assembly'],[[r.run,r.model,r.assembly] for r in Q.itertuples()],[112,208,190])
p('All seven run accessions, original download URLs and study identifiers are in config/reads.csv. Two NextSeq runs share PRJEB85433 and a release date; they are candidate batch mates, not a proven shared flow-cell lane. The NovaSeq phytoplasma-associated run is registered under its plant host, illustrating why organism labels alone do not reliably join read and assembly records.')
p('The input archive files were retrieved on 29 September 2026 and reused locally; this new analysis and BLAST verification were run on 30 September. The fixed manifests and input checksums record what was used. scripts/fetch.py retrieves raw data; scripts/metadata.py re-fetches versioned metadata into a separate directory so a live refresh cannot silently replace the study snapshot.',True)
page()

title('3. Methods: a short and visible decision rule')
sub('Assembly screen')
p('Every contig is scanned in both orientations for four adapter motifs (TruSeq R1/R2, Nextera and BGI). Exact 19-21 bp matches avoid an arbitrary approximate-alignment score, but miss mutated and truncated adapters. All contigs of at least 1 kb are also mapped with minimap2 through mappy, using the asm5 preset, to the foreign panel and UniVec_Core. A retained match has at least 200 matching bases, identity at least 95% and MAPQ at least 20. Shorter contigs still receive the adapter scan.')
p('Minimap2 finds shared minimizer seeds, chains them and extends promising chains into alignments. It is a heuristic, so its speed does not imply perfect sensitivity. Matches to the same reference are combined by the union of query intervals; this prevents a rotated circular mitochondrial sequence from being missed when it aligns in two pieces.')
table(['Evidence','Score'],[['At least 80% of the contig covered by the foreign reference','+2'],['Alignment identity at least 98%','+1'],['Contig GC differs from assembly-wide GC by at least 12 percentage points','+1'],['Foreign match has at least 500 bp of sequence on both sides','-3']],[450,60])
p('A score of at least 3 is a suspected foreign contaminant. A foreign block with long flanks is retained as an HGT-or-construct candidate; other matches remain uncertain. The score is a prioritization rule, not a trained probability. Vector matches are listed separately because vectors contain natural bacterial segments and engineered genomes may intentionally contain them.')
sub('Independent database check and baseline')
p('The strongest foreign candidates are submitted by script to NCBI BLAST using megablast, expectation threshold 1e-10 and up to 20 hits. Self-accession hits are excluded. BLAST uses word seeds and extension to score local similarity; identity and query coverage are reported together. A 1 kb BLAST hit identifies that window, while whole-contig support comes from the panel alignment.')
p('On 100 real 1 kb native windows, the adapter screen returned zero flags. Adding a perfect Nextera motif was detected in 100/100 cases. One substitution in that motif reduced detection to 0/100. This measured weakness is why the observed archive rate must not be presented as total contamination prevalence.')
page()

title('4. Archive results and platform / era patterns')
figure('archive.png','Figure 1. Assemblies with each screening signal and the main flag rate by declared platform. Categories overlap. Vector and HGT/construct matches (*) are excluded from the main flag rate.',200)
figure('years.png','Figure 2. Main flag rate by assembly release year. Small early-year denominators and project composition make year-to-year comparisons unstable.',150)
p(f'The main rate is {H["flagged_assemblies"]}/{H["assemblies"]} = {pct(H["flag_rate"])}. The Wilson interval is {pct(H["flag_rate_ci95"][0])} to {pct(H["flag_rate_ci95"][1])}. This interval summarizes binomial count uncertainty; correlated samples and the fixed archive census mean it is not a correction for ascertainment bias or unknown false negatives.')
rows=[]
for _,r in rates[(rates.grouping=='era')&(rates.category=='flagged')].iterrows():rows.append([r['group'],int(r.n),int(r.flagged),pct(r.rate)])
table(['Release era','Assemblies','Flagged','Rate'],rows,[180,110,110,110])
p('Full category-by-year, taxon and platform counts and denominators are supplied in results/rates.csv. The platform comparison is observational: it mixes species, laboratories, MAGs and cultured isolates. More flags in short-read assemblies do not establish that an instrument alone caused contamination. Adapter motifs reflect library chemistry; assembly fragmentation also changes how many exposed contig ends can retain them.',True)
page()

title('5. Short reads: quality and technical sequence')
p('Four Illumina runs contribute the first 10,000 pairs each (20,000 individual reads per run). The analysis records quality by cycle, GC distribution, exact adapter matches, duplicate sequences, poly-G tails and trinucleotide entropy. A low-complexity flag is entropy below 3 bits; it is descriptive because AT-rich bacterial reads can be naturally simple. A run-specific GC histogram is used instead of calling every second peak a contaminant.')
figure('short_reads.png','Figure 3. Raw base quality across sequencing cycles and read GC distributions for the MiSeq and host-associated NovaSeq examples. R1 and R2 contribute to the same cycle summary.',200)
short=Q[Q.platform=='Illumina']
table(['Run','Adapter reads','Trimmed','Discarded','Mean Q before / after'],[[r.run,pct(r.adapter_fraction),pct(r.trimmed_fraction),pct(r.discarded_fraction),f'{r.mean_phred:.1f} / {r.mean_retained_phred:.1f}'] for r in short.itertuples()],[112,92,85,85,136])
p('Before mapping, the first exact adapter is removed, 3-prime bases below Phred 20 are trimmed, and reads shorter than 50 bp are discarded. Q20 corresponds to a nominal 1% base error probability; the 50 bp minimum retains enough sequence for reasonably specific placement. The exact adapter rule remains conservative and misses partial adapters. Reported raw-Q averages are read-weighted and retained-Q averages are base-weighted, so their difference is descriptive rather than a controlled effect estimate.')
r=Q[Q.run=='SRR9681681'].iloc[0]
p(f'The MiSeq example SRR9681681 has adapter sequence in {pct(r.adapter_fraction)} of reads. Its adapter burden is much larger than the other selected Illumina runs. The host-associated NovaSeq example maps poorly to the bacterial assembly and our small panel. Plant nuclear sequence is a plausible explanation, but unmapped reads cannot be assigned a taxon without a broader search.')
p('Duplication is measured by repeated sequences, not by optical duplicate coordinates. Poly-G means at least ten terminal G bases; this is a warning for two-colour chemistry, not proof of an instrument error. No read-level rates are extrapolated to all 3,006 assemblies.',True)
page()

title('6. Long reads and the index-hopping check')
figure('platforms.png','Figure 4. Read length and observed alignment difference. Different organisms and sequencing generations are mixed; this is a descriptive comparison, not a paired technology benchmark.',200)
long=Q[Q.platform!='Illumina']
table(['Run / technology','Median bp','Identity','Native + foreign splits'],[[r.run+' / '+r.platform,f'{r.median_length:,.0f}',pct(r.median_identity),f'{r.foreign_native_splits}/2000'] for r in long.itertuples()],[205,85,85,135])
p('Each long-read run contributes the first 2,000 reads. ONT uses map-ont and the PacBio Sequel IIe run uses map-hifi. Error rates are measured against the selected reference, so they include strain differences and alignment effects as well as sequencing errors. Very high archived PacBio quality values must not be compared directly with old ONT qualities as if all basecallers were calibrated alike.')
p('The chimera check requires two alignments with MAPQ at least 20, at least 300 query bases each and at most 50 bp overlap; one must map to the native assembly and the other to the foreign panel. No such read was detected in these three prefixes. With 0/2,000 per run, the Wilson upper bound is about 0.19%, conditional on this detector and sampling scheme. Same-organism rearrangements, palindromes and contaminants absent from the panel are not measured.')
sub('Candidate pooled libraries')
p('In each of ERR14335635 and ERR14335636, 10,000 pairs were mapped competitively to both assemblies. Both mates must uniquely favor the other sample at MAPQ at least 30, identity at least 98% and query coverage at least 90%. The result was zero cross-assigned pairs in both runs (upper Wilson bound about 0.038% per run).')
p('This is a screen for cross-sample sequence, not a direct measurement of index hopping. Index reads, lane assignments and the original sample sheet were unavailable. Similar biological sequence, laboratory mixing and genuine index hopping cannot be separated from demultiplexed FASTQ alone. The two negative results therefore do not establish that index hopping was absent.')
page()

title('7. Error analysis and case review')
p('A seeded random sample selected one representative hit per assembly: 20 adapter assemblies and 10 vector assemblies, plus all three foreign candidates. The saved sample was examined case by case using the surrounding sequence, contig position, BioSample descriptions and non-self BLAST results. Every decision and its reason is recorded in results/review_decisions.csv. This is documented evidence adjudication, not a blinded human-validation study or laboratory truth set.')
review=json.loads((R/'review_summary.json').read_text())
table(['Screen','Reviewed','Supported','False label','Unresolved'],[[c,v['reviewed'],v['supported'],v['false'],v['unresolved']] for c,v in review.items()],[150,90,90,90,90])
p('Nineteen adapter cases had strong contextual support, including extended adapter sequence, an exposed contig end or an adapter immediately before an assembly gap. One isolated internal Nextera match remained unresolved. Among the 19 resolved cases, the observed false-discovery fraction is 0/19, with a 95% Wilson upper bound of 16.8%. Across all 20 inspected cases, unresolved-case bounds are 0-5%; these bounds are not a confidence interval for the whole archive.')
p('All three foreign candidates were supported by both whole-contig alignment and an independent non-self BLAST window. Zero rejected cases out of three is weak error-rate evidence: the Wilson upper bound is 56.2%. This small denominator prevents a claim that the foreign detector has near-perfect specificity.')
p('Six of ten reviewed vector cases were GM12 strains whose BioSample descriptions explicitly stated construction or deletion-mutant work. Their internal E. coli lac-operon segments are compatible with intentional constructs. Treating these as accidental contamination would give at least 6/10 false labels in that sample. Four other partial vector matches remained unresolved because shuttle vectors can carry natural microbial sequence. This is the practical reason vector matches are excluded from the main count.')
sub('Two different error rates')
p('The handbook calls for a false-positive estimate from flagged records. Strictly, false labels divided by resolved flagged records estimate a false-discovery fraction; a conventional false-positive rate requires a representative set of truly clean negatives. The 100 native-window check is a small detector baseline and is not such a genome-wide truth set.')
sub('Sensitivity to the score')
v=pd.read_csv(R/'alignment_hits.csv');f=v[v.category=='Foreign'];p('Changing the foreign-call score threshold from 3 to 2 or 4 yields '+', '.join(f'{t}: {f[f.score>=t].accession.nunique()} assemblies' for t in [2,3,4])+'. These are sensitivity checks on already retained alignments, not a re-benchmark of the identity or MAPQ filters.')
page()

title('8. Ten findings for a database curator')
p('Each row is a review request tied to a versioned record. Coordinates are 1-based and inclusive. All matches are available in adapter_hits.csv and alignment_hits.csv; the three foreign cases also have query sequences, database responses and nt_hits.csv.')
top=pd.read_csv(R/'curator_top10.csv')
for i,r in enumerate(top.itertuples(),1):
 sub(f'{i}. {r.accession} - {r.category}')
 p(f'<b>{r.contig}:{r.start}-{r.end}</b><br/>{escape(r.evidence)}. {escape(r.action)}',True)
p('The ten candidates cover three foreign-DNA cases and seven adapter-rich assemblies. They are prioritized for clarity and evidence, not chosen to imply that every possible category produced a confirmed mistake. There was no demonstrated sample swap or index-hopping event in this analysis.',True)
page()

title('9. Reproducibility, limits and conclusion')
p('Install the pinned environment and run <b>snakemake --cores 2</b> for the bundled real mini-dataset. Run <b>snakemake --cores 2 --config mode=full</b> to fetch the fixed archive inputs and regenerate the headline result. BLAST responses are keyed by the exact query SHA-256 and can be replayed offline. scripts/check.py verifies the mini-dataset counts and the adapter-spike baseline.')
scan=H['scan_runtime']['seconds'];align=H['align_runtime']['seconds'];reads=H['reads_runtime']['seconds'];rss=H['align_runtime']['peak_rss_os_units']/1024**2
p(f'Measured component times on this macOS machine were {scan:.1f} s for scanning, {align:.1f} s for panel/vector alignment and {reads:.1f} s for read QC. Alignment peak process memory was {rss:.0f} MiB (macOS reports ru_maxrss in bytes). Download time, environment setup and the NCBI queue are excluded. FASTA input totals 3.67 Gb; the intermediate query file is about 3.4 GiB. Genomes are streamed one assembly at a time; the candidate metadata table is kept in memory.')
p('At the same observed throughput, ten times as much sequence would take roughly ten times the scan/alignment time, assuming similar repeats and reference size. This is an extrapolation, not a measured tenfold run. More cores can overlap independent workflow steps; the current per-contig loop is serial. Cost is mainly download, CPU time and storage. Actual monetary cost cannot be reconstructed from archive records; platform price claims would be speculative.')
sub('Limits that change the interpretation')
p('The foreign panel has only four references; nuclear host DNA, divergent organisms and short vector fragments can be missed. Exact motifs miss damaged adapters. GC can vary naturally. HGT and constructs cannot be proved from contig context alone. All raw reads are prefixes, not random samples. Species, laboratories and eras confound the platform comparison. The direct metadata join misses parent-sample links. core_nt omits some large eukaryotic chromosomes. Case review was small and not an independent blinded validation. The measured 8.18% flag rate is therefore a narrow triage result.')
sub('Conclusion and technical responsibilities')
p('Simple screening found useful, traceable candidates at archive scale. The strongest foreign cases were E. coli and mouse mitochondrial DNA, while technical adapter residue was more frequent. A curator should inspect sequence context and raw reads before editing a record. Abylaikhan Torekhan and Kuanysh Bakizhan both worked on the code and checked the results. Both partners are responsible for the final report and for explaining the full analysis at the defence.')
sub('References')
refs=[('NCBI VecScreen and interpretation','https://www.ncbi.nlm.nih.gov/tools/vecscreen/interpretation/'),('NCBI nucleotide database scope','https://blast.ncbi.nlm.nih.gov/doc/blast-help/blastdatabases.html'),('Illumina adapter sequences','https://support-docs.illumina.com/SHARE/AdapterSequences/Content/Nextera_Illumina-Sequences.htm'),('Salter et al. 2014. Reagent and laboratory contamination','https://doi.org/10.1186/s12915-014-0087-z'),('Li 2018. Minimap2','https://doi.org/10.1093/bioinformatics/bty191'),('Illumina: index hopping','https://supportassets.illumina.com/techniques/sequencing/ngs-library-prep/multiplexing/index-hopping.html')]
for name,url in refs:p(f'<link href="{url}" color="#20639b">{name}</link>',True)

def footer(canvas,doc):
    canvas.setFont('Helvetica',8);canvas.setFillColor(colors.HexColor('#637589'))
    canvas.drawString(42,24,'CSI Genome | Torekhan & Bakizhan | BDA-2405');canvas.drawRightString(A4[0]-42,24,str(doc.page))
SimpleDocTemplate(str(ROOT/'report/CSI_Genome_Report_draft.pdf'),pagesize=A4,rightMargin=42,leftMargin=42,topMargin=36,bottomMargin=38,
  title='CSI Genome - Project 07',author='Abylaikhan Torekhan; Kuanysh Bakizhan').build(story,onFirstPage=footer,onLaterPages=footer)
print('Report saved')
