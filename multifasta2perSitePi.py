#!/usr/bin/env python3
# coding: utf-8
# Author : Diana Aguilar


import argparse
import pandas as pd
import Bio
from Bio import SeqIO
import itertools
import matplotlib.pyplot as plt

#Read file
parser = argparse.ArgumentParser()
parser.add_argument("fastaname", type=str, help="multialignment fasta")
args = parser.parse_args()
seq_var=list(Bio.SeqIO.parse(args.fastaname, "fasta"))

#Number of comparisons
num=int((len(seq_var)*(len(seq_var)-1))/2)
print(num," pairs of sequences")

#Do per base comparison
pairwise=[]
for seq1,seq2 in list(itertools.combinations(seq_var, 2)):
    compare=[int(x!=y) for x,y in zip(seq1.seq,seq2.seq)]
    pairwise.append(compare)
    
#Average per base pair and save to file
pi_persite_df=pd.DataFrame(pairwise).mean(axis=0)
pi_persite_df.to_csv(args.fastaname+".Pairwise",sep="\t",index=False,header=False)

#Plot
pi_persite=list(pi_persite_df)
plt.figure(figsize=(10, 5), dpi=80)
plt.scatter(range(0,len(pi_persite)),pi_persite)
plt.xlabel("Position")
plt.ylabel("pi")
plt.savefig(args.fastaname+"image.png")
