from Bio import SeqIO
import math
import argparse

def Clean(fasta, outfile, min_len=None, max_len=None):
    min_len=min_len if min_len is not None else 0
    max_len=max_len if max_len is not None else math.inf
    seqs={}
    for rec in SeqIO.parse('pfam_unique.fasta',format='fasta'):
        length=len(rec.seq)
        if 'X' in rec.seq or 'B' in rec.seq or (length<min_len or length>max_len):
            continue
        # id ='|'.join(rec.id.split('|')[:2])
        id=rec.id
        seqs[id]=str(rec.seq)

    # SeqIO.write(seqs, 'pfam_unique_filtered.fasta', format='fasta')
    with open('pfam_unique_filtered.fasta', 'w') as f:
        for id, seq in seqs.items():
            print(f'>{id}', file=f)
            print(seq, file=f)

if __name__=='__main__':
    parser=argparse.ArgumentParser(description="Basic filtering of fasta file by length, ambiguous characters. Also trims header to the ID.")
    parser.add_argument('fasta', help='input fasta file')
    parser.add_argument('output', help='output fasta file')
    parser.add_argument('-min_len', type=int, help='Exclude sequences less than this length')
    parser.add_argument('-max_len', type=int, help='Exclude sequences grater than this length')
    args=parser.parse_args()

    Clean(args.fasta, args.output, args.min_len, args.max_len)
