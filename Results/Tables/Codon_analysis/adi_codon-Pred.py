import pandas as pd
import csv

prefix = 'adi_'
Conds = ['D0','D4','D8'] 
reps = ['a','b','c']

#data file
folder = "/Data_2/Jun/Adipocytes/rpf/novaseq/"
SPs = [Cond+rep for Cond in Conds for rep in reps]
filetype = ".rep.codons.txt"

#%%
def codon2index(codon):
    if 'N' in codon : return -1
    a = {'U':0, 'C':16, 'A':32, 'G':48}
    b = {'U':0, 'C':4,  'A':8,  'G':12}
    c = {'U':0, 'C':1,  'A':2,  'G':3}
    return a[codon[0]]+b[codon[1]]+c[codon[2]]

#%%
Pout_w = pd.ExcelWriter(f'{prefix}P-codon-Pred.xlsx')
Aout_w = pd.ExcelWriter(f'{prefix}A-codon-Pred.xlsx')

for SP in SPs:
    P_result = {}
    A_result = {}
    P_result['transcript_id'] = [a+b+c for a in ['U','C','A','G'] for b in ['U','C','A','G'] for c in ['U','C','A','G']]
    A_result['transcript_id'] = [a+b+c for a in ['U','C','A','G'] for b in ['U','C','A','G'] for c in ['U','C','A','G']]
    
    with open(folder+SP+filetype) as f:
        tr = csv.reader(f, delimiter='\t')
        next(tr, None)
        for row in tr:
            # transcript_id, psite asite, codon-psite, codon-asite, aa-psite, aa-asite, length, frame, reltostart-asite, reltostop-asite
            if P_result.get(row[0])==None:
                P_result[row[0]] = [0]*64
                A_result[row[0]] = [0]*64
                
            if codon2index(row[3]) != -1 : P_result[row[0]][codon2index(row[3])] += 1
            if codon2index(row[4]) != -1 : A_result[row[0]][codon2index(row[4])] += 1

    print(f'{SP} done')

    P_result = pd.DataFrame(P_result).T
    A_result = pd.DataFrame(A_result).T
    
    P_result.to_excel(Pout_w, sheet_name=SP, header=False)
    A_result.to_excel(Aout_w, sheet_name=SP, header=False)

Pout_w.close()
Aout_w.close()
