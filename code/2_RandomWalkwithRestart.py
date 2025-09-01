#/usr/bin/python3
import numpy as np
import pandas as pd
import time
import random
from scipy.sparse import csr_matrix
from scipy.sparse import csc_matrix
from scipy import sparse
import multiprocessing as mp
import sys
import re
import torch
import os
os.chdir("/XXX")
os.getcwd()
def my_rand_walk(p0,W,r):
    pt = p0
    delta = 1
    while(delta > 1e-10):
        pt1 = (1-r) * np.dot(W, pt) + r * p0
        delta = sum(abs(pt1-pt))
        pt = pt1
    return(pt)

backgroundnetwork_file = "network.txt"
r = 0.7
p_matrix_out = "".join(("/XXX/", str(r), "/"))
net = pd.read_csv(backgroundnetwork_file, sep = '\t', header = None, names = ['n1','n2'], dtype = {'n1':str, 'n2':str})
nodes = list(np.loadtxt("nodes.txt", dtype = str))
m = np.zeros((len(nodes),len(nodes)))
for i in range(0, net.shape[0]):
    n1 = nodes.index(net.iloc[i, 0])
    n2 = nodes.index(net.iloc[i, 1])
    m[n1,n2] = 1
    m[n2,n1] = 1

W = m/sum(m)
p_matrix = np.eye(len(nodes), len(nodes))## diagonal matrix
mRNA_index = list(pd.read_csv("mRNA_index.txt", sep = "\t", header = None)[0])
miRNA_index = list(pd.read_csv("miRNA_index.txt", sep = "\t", header = None)[0])
lncRNA_index = list(pd.read_csv("lncRNA_index.txt", sep = "\t", header = None)[0])
circRNA_index = list(pd.read_csv("circRNA_index.txt", sep = "\t", header = None)[0])
ncRNA_index = miRNA_index + lncRNA_index + circRNA_index
mRNA = list(pd.read_csv("mRNA.txt", sep = "\t", header = None, dtype = str)[0])
miRNA = list(pd.read_csv("miRNA.txt", sep = "\t", header = None)[0])
lncRNA = list(pd.read_csv("lncRNA.txt", sep = "\t", header = None)[0])
circRNA = list(pd.read_csv("circRNA.txt", sep = "\t", header = None)[0])
ncRNA = miRNA + lncRNA + circRNA
p_matrix = p_matrix[:,ncRNA_index]
for i in range(0,127):
    print(i)
    out = "".join(("out",str(i)))
    final_matrix = np.apply_along_axis(my_rand_walk,0,p_matrix[:,(100*i):(100*(i+1))],W,r)
    np.savetxt(out, final_matrix, fmt='%.18e', delimiter='\t', newline='\n', header='', footer='', encoding=None)

i = 0
file = "".join(("out",str(i)))
RWR_mat = np.loadtxt(file)
for i in range(1,127):
    file = "".join(("out",str(i)))
    infile = np.loadtxt(file)
    RWR_mat = np.append(RWR_mat,infile,axis = 1)

np.savetxt("out", RWR_mat, fmt='%.18e', delimiter='\t', newline='\n', header='', footer='', encoding=None)
RWR_mat = pd.read_csv("out", sep = '\t', header = None)

for i in range(0,RWR_mat.shape[1]):
    outfile = pd.DataFrame({"index": mRNA, "score": RWR_mat.iloc[mRNA_index, i]})
    np.savetxt("_".join((ncRNA[i],"finalScore.rnk")), outfile, fmt='%s\t%.18e', newline='\n', header='', footer='', encoding=None)

