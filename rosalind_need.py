from Bio import Entrez
import numpy as  np
import pandas as pd
Entrez.email = 'mail@mail.com'

def importer(ids : str):
    ids = ids.split()
    seqs = []
    for i in ids:
        print(i+' import started')
        handle = Entrez.efetch(db='nucleotide',id=i, rettype = 'fasta', retmode= 'text')
        record = handle.read().split('\n',1)
        handle.close()
        print(i+' import finished')
        record[1]=record[1].replace('\n','')
        seqs.append(record[1])
    return seqs



table = pd.read_csv('dnafull.csv')
table = table.set_index('O')


def needleman_wunch(seqs]):
    '''
    Needleman-Wunsch algorithm based on dnafull table,
    gap open penalty is 10, gap extension penalty is 1 (affine gap penalty)

    Parameters
    ----------
    seqs : an list of 2 seqs for pairwise alignment

    Returns
    -------
    final_array : np.array with values of alignment
    [-1][-1] index denotes final score of alignment

    '''
    seq1,seq2 = seqs
    gap_open = -9
    gap_ext = -1
    len_seq1 = len(seq1)+1
    len_seq2 = len(seq2)+1
    
    G = np.zeros((len_seq1,len_seq2))
    F = np.zeros(shape=(len_seq1,len_seq2))   #сначала вертикаль потом горизонталь
    E = np.zeros(shape=(len_seq1,len_seq2))
    V = np.zeros(shape=(len_seq1,len_seq2))
    G[0,:] = 0
    G[:,0] = 0
    G[0,0] = -1*np.inf
    E[:,0] = [0]+ [gap_open+i*gap_ext for i in range(1,len_seq1)]
    E[0,:] = 0
    E[0,0] = -1*np.inf
    F[0,:] = [0]+ [gap_open+i*gap_ext for i in range(1,len_seq2)]
    F[:,0] = 0
    F[0,0] = -1*np.inf
    V[0,:] = [0]+ [gap_open+i*gap_ext for i in range(1,len_seq2)]
    V[:,0] = [0]+ [gap_open+i*gap_ext for i in range(1,len_seq1)]
    V[0,0] = 0
    for i in range(1,len_seq1):
        for j in range(1,len_seq2):
            E[i,j] = max(E[i,j-1]+gap_ext, V[i,j-1]+gap_open+gap_ext)
            F[i,j] = max(F[i-1,j]+gap_ext, V[i-1,j]+gap_open+gap_ext)
            G[i,j] = V[i-1,j-1]+ table[seq1[i-1]][seq2[j-1]]
            V[i,j] = max(G[i,j], E[i,j], F[i,j])
    final_array=V
    return final_array


z = needleman_wunch(importer('JX205496.1 JX469991.1'))

print(z[-1][-1])
