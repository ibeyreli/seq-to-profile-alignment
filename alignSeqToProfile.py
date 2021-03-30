# -*- coding: utf-8 -*-
"""
CS481 Fall 2019 Homework 4
Sequence to Profile Alignment
Source Code File
Created on Thu Nov 7 18:09:11 2019
Revised on Thu Nov 26 19:50:39 2019
@author: İlayda Beyreli 201801130
"""
import sys
import time
from itertools import chain

def readfa_multi(file_name):
    # Reading .FA files into lists
    seqs  = []
    names = []
    with open(file_name, "r") as fin:
        data =fin.read().splitlines()
    if ">" in data[0]:
        names = data[0][2:]+" "
        seqs = list(chain.from_iterable(data[1]))
    else:
        while data:
            seq1 = data[0]
            s = seq1.index(" ")
            n1 = seq1[:s]
            seq1 = seq1[s+1:]
            seq1 = list(chain.from_iterable(seq1))
            seqs.append(seq1)
            names.append(n1)
            data = data[1:]
        seqs = list(zip(*seqs))
        seqs = [list(seq) for seq in seqs]
    return names, seqs

def writefa_multi(s,names,seqs,aseqp,name,aseq,file_name):
    seqs = list(zip(*seqs))
    seqs = [list(seq) for seq in seqs]
    s = 0
    s = [s for s in range(len(aseqp)) if aseqp[s] == "-"]
    with open(file_name, "w+") as fout:
        i=0
        while names:
            for j in s:
                seqs[i].insert(j, "-")
            seqs[i] = "".join(seqs[i])

            fout.write(names.pop(0)+" "+seqs[i]+"\n")
            i+=1
        fout.write(name+" "+aseq)
    return True

def zeros(r,c):
    # Matrix as list of lists 
    return [[0 for i in range(c)] for i in range(r)]

def matprint(matrix):
    for i in range(len(matrix)):
            print(matrix[i],"\n")
    return 0

def profiler(seqs,alphabet):
    n = len(seqs)
    profile = zeros(n,5)
    for i in range(n):
        for j in range(5):
            profile[i][j] = seqs[i].count(alphabet[j])/5.0
    return profile

def naive_score(c1,c2,p,match,gapop,miss):
    # Prob is the probabilities for that significant location in profile
    # Scoring Function, given in the assignment
    if c2 =='-' or c1 =='-':
        return gapop*p
    elif c2 == c1:
        return match*p
    else:
        return miss*p #missmatch
    
def profile_score(prf,c1,alp,match,gapop,miss):
    score = 0
    for i in range(len(alp)):
        c2 = alp[i]
        p = prf[i]
        score += naive_score(c1,c2,p,match,gapop,miss)
    return score
    
def argmax(a1,a2,a3):
    switcher={0:'d',
            1:'u',
            2:'l'}
    l = [a1,a2,a3]
    m = max(a1,a2,a3)
    m = l.index(m)
    return switcher.get(m)
    
def naive_prof_alignment(seq,profile,match,gapop,miss):
    # Needleman-Wunsh & Smith-Whitterman Naive Alignment Algorithms
    alp = ['A','T','G','C','-']
    # Initialize table
    n = len(profile)#i columns
    m = len(seq)# j rows
    S = zeros(m+1,n+1)
    L = zeros(m+1,n+1)

    # Initialize according to global alignment
    for i in range(1,n+1):
        S[0][i] += S[0][i-1]+profile_score(profile[i-1],'-',alp,match,gapop,miss)
    for j in range(1,m+1):
        S[j][0] += S[j-1][0]+naive_score(seq[j-1],'-',1,match,gapop,miss)
    # Fill the scoring table
    #matprint(S[0])
    for j in range(1,m+1):
        for i in range(1,n+1):
            match = S[j-1][i-1]+profile_score(profile[i-1],seq[j-1],alp,match,gapop,miss)
            del1 = S[j-1][i]+naive_score(seq[j-1],'-',1,match,gapop,miss) 
            in1 = S[j][i-1]+profile_score(profile[i-1],'-',alp,match,gapop,miss)
            S[j][i] = max(match,del1,in1)
            L[j][i] = argmax(match,del1,in1)
    # Score for global alignment
    alignscore = S[-1][-1]
    # matprint(S) # Debugging
    # matprint(L) # Debugging
    # Traversing
    aseqp=[]
    aseq=[]
    i = n# [j for j in range(n+1) if alignscore in S[j]][0]
    j = m #S[j].index(alignscore) 
    while j > 0 and i > 0:
        if S[j][i] == 0 and alignscore !=0:
            break
        # print(L[i][j]) # Debugging
        if L[j][i] == 'd':
            aseqp.append("*")
            aseq.append(seq.pop(-1))
            i-=1
            j-=1
        elif L[j][i] == 'u':
            aseqp.append('-')
            aseq.append(seq.pop(-1))
            j-=1
        elif L[j][i] == 'l':
            aseqp.append("*")
            aseq.append('-')    
            i-=1
        else:
            i-=1
            j-=1
    # print(aseq1,aseq2) # Debugging
    aseqp.reverse()
    aseqp=''.join(aseqp)
    aseq.reverse()    
    aseq=''.join(aseq)
    return alignscore, aseqp,aseq
"""
# Test & Debug
# aligned_sequences.aln  # seq.fa
match = 1
gapop = -2
miss = -1
alphabet = ['A','T','G','C','-']
names, seqs = readfa_multi("aligned_sequences.aln")
name, seq = readfa_multi("seq.fa")
print(''.join(seq))
prf = profiler(seqs,alphabet)
#print(names, seqs)
print(prf)
s,aseqp,aseq= naive_prof_alignment(seq,prf,match,gapop,miss)
print(s,aseq,aseqp)
print(len(aseq))
# sequence  CTAGA---TAATTGGAGATGATCAAAT-TTATAT
"""
# Main Function to Run From Terminal
# alignSeqToProfile --fasta seq.fasta --aln aligned_sequences.aln --out seq.aln
# --gap ${gap_penalty} --match ${match_score} --mismatch ${mismatch_penalty}
try:
    e = 13
    if len(sys.argv) < e:
        raise IndexError
    if sys.argv[0] != "alignSeqToProfile.py":
        arg = sys.argv[1]
        raise NameError
    seq = [sys.argv[i] for i in range(1,len(sys.argv)) if sys.argv[i-1]=="--fasta"][0]
    aln = [sys.argv[i] for i in range(1,len(sys.argv)) if sys.argv[i-1]=="--aln"][0]
    out = [sys.argv[i] for i in range(1,len(sys.argv)) if sys.argv[i-1]=="--out"][0]
    gapop = [int(sys.argv[i]) for i in range(1,len(sys.argv)) if sys.argv[i-1]=="--gap"][0]
    match = [int(sys.argv[i]) for i in range(1,len(sys.argv)) if sys.argv[i-1]=="--match"][0]
    mis = [int(sys.argv[i]) for i in range(1,len(sys.argv)) if sys.argv[i-1]=="--mismatch"][0]
except IndexError:    
    print("Insufficient number of arguments! Expected:",e, " Passed:",len(sys.argv)-1)
except NameError:
    print("Unknown Argument", arg)
    
alphabet = ['A','T','G','C','-']
names, seqs = readfa_multi(aln)
name, seq = readfa_multi(seq)    
prf = profiler(seqs,alphabet)
s,aseqp,aseq= naive_prof_alignment(seq,prf,match,gapop,mis)
Done = writefa_multi(s,names,seqs,aseqp,name,aseq,out)
