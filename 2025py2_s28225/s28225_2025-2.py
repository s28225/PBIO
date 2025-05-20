from Bio import Entrez, SeqIO;import csv,matplotlib.pyplot as p
class R:
 def __init__(s,e,k):Entrez.email,Entrez.api_key=e,k
 def search(s, t, min=None, max=None):s.w,s.q=Entrez.read(Entrez.esearch(db="nucleotide",term=f"txid{t}[Organism]{f' AND {min}:{max}[SLEN]' if min and max else f' AND {min}:[SLEN]' if min else f' AND 0:{max}[SLEN]' if max else ''}",usehistory="y"))["WebEnv"], Entrez.read(Entrez.esearch(db="nucleotide",term=f"txid{t}[Organism]{f' AND {min}:{max}[SLEN]' if min and max else f' AND {min}:[SLEN]' if min else f' AND 0:{max}[SLEN]' if max else ''}",usehistory="y"))["QueryKey"]
 def fetch(s,n=10):return list(SeqIO.parse(Entrez.efetch(db="nucleotide",rettype="gb",retmode="text",retstart=0,retmax=n,webenv=s.w,query_key=s.q),"genbank"))
 def to_csv(s,r,f):csv.writer(open(f,"w",newline='',encoding='utf-8')).writerows([["Accession","Length","Description"]]+[[x.id,len(x.seq),x.description]for x in r])
 def plot(s,r,f):r.sort(key=lambda x:-len(x.seq));p.plot([x.id for x in r], [len(x.seq) for x in r], 'o-b');p.xticks(rotation=45, ha='right');p.tight_layout();p.grid();p.savefig(f);p.close()
def main():
 e,k,t=input("Email:"),input("Api key:"),input("TaxID:")
 a,b=[int(x)if x.isdigit()else None for x in[input("Min:"),input("Max:")]]
 r=R(e,k);r.search(t,a,b);d=r.fetch()
 SeqIO.write(d,f"taxid_{t}_sample.gb","genbank");r.to_csv(d,f"taxid_{t}_summary.csv");r.plot(d,f"taxid_{t}_plot.png")
if __name__ == "__main__":main()