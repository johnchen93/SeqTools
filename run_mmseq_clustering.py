import argparse
import subprocess
import tempfile
import os
import pandas as pd
from Bio import SeqIO

def Cluster(fasta_files, id, label):
    with tempfile.tempdir() as tmpdir:
        db_path=os.path.join(tmpdir, 'db')
        clu_path=os.path.join(tmpdir, 'clu')
        clu_tmp=os.path.join(tmpdir,'clu_tmp')
        clusters_out=f'{label}_{id}.clust'
        c=id/100
        cluster_param=['--min-seq-id',f'{c:.2f}']
        fasta_out=f'{label}_{id}.fasta'
        subprocess.call(['mmseqs','createdb']+fasta_files+[db_path])
        subprocess.call(['mmesqs','linclust',db_path, clu_path, clu_tmp]+cluster_param)
        subprocess.call(['mmseqs','createtsv',db_path, db_path, clu_path, clusters_out])

        clusters=pd.read_csv(clusters_out, sep='\t', header=None, names=['rep','member'])
        reps=set(clusters['rep'].unique())

        seqs=[] 
        for file in fasta_files:
            for fasta in SeqIO.parse(file, format='fasta'):
                if fasta.id in reps:
                    seqs.append(fasta)
        SeqIO.write(fasta_out, format='fasta')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('fastas', nargs='+', help='list of fastas to process')
    parser.add_argument('cutoff', type=int, help='cutoff to cluster at where 100 is 100%')
    parser.add_argument('label', help='Label for output, extra suffixes (including id cutoff) will be added.')
    args=parser.parse_args()

    Cluster(args.fastas, args.cutoff, args.label)