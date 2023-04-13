import argparse
import subprocess
import tempfile
import os
import pandas as pd
from Bio import SeqIO

def Cluster(fasta_files, id, clusters_out, full_clust=False):
    with tempfile.TemporaryDirectory() as tmpdir:
        db_path=os.path.join(tmpdir, 'db')
        clu_path=os.path.join(tmpdir, 'clu')
        clu_tmp=os.path.join(tmpdir,'clu_tmp')
        func='cluster' if full_clust else 'linclust'

        c=id/100
        cluster_param=['--min-seq-id',f'{c:.2f}']
        subprocess.call(['mmseqs','createdb']+fasta_files+[db_path])
        subprocess.call(['mmseqs',func, db_path, clu_path, clu_tmp]+cluster_param)
        subprocess.call(['mmseqs','createtsv',db_path, db_path, clu_path, clusters_out, '--full-header'])

def WriteCluster(fasta_files, clusters_out, fasta_out):
    clusters=pd.read_csv(clusters_out, sep='\t', header=None, names=['rep','member'])
    reps=set(clusters['rep'].unique())

    with open(fasta_out, 'w') as f:
        for file in fasta_files:
            for fasta in SeqIO.parse(file, format='fasta'):
                # print(fasta.id)
                if fasta.id in reps:
                    # print(fasta.id)
                    f.write(f'>{fasta.id}\n{fasta.seq}\n')

if __name__=='__main__':
    parser=argparse.ArgumentParser()
    parser.add_argument('fastas', nargs='+', help='list of fastas to process')
    parser.add_argument('cutoff', type=int, help='cutoff to cluster at, where 100 is 100%, 95 is 95%, etc.')
    parser.add_argument('label', help='Label for output, extra suffixes (including id cutoff) will be added.')
    parser.add_argument('-force', action='store_true', help='Flag. Force redo even if file already exists')
    parser.add_argument('-no_fasta', action='store_true', help='Flag. If set, do not return fasta file')
    parser.add_argument('-full_clust', action='store_true', help='Flag. If set, use cluster instead of linclust')
    args=parser.parse_args()

    clusters_out=f'{args.label}_{args.cutoff}.clust'
    if args.force or not os.path.exists(clusters_out):
        Cluster(args.fastas, args.cutoff, clusters_out, args.full_clust)
    if not args.no_fasta:
        fasta_out=f'{args.label}_{args.cutoff}.fasta'
        WriteCluster(args.fastas, clusters_out, fasta_out)