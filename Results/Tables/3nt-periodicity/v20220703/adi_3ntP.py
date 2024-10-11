import pandas as pd
import pysam

col_list = ["Transcript_ID","Gene_ID","Gene_name","5UTR_length","CDS_length","3UTR_length"]
P_gtf = pd.read_csv('/Data_1/Daehwa/python/GTF_Parsing/v2.2/gencode.vM27.annotation.gtf/v20220414/Processed_gtf.tsv', sep='\t', usecols=col_list)
repre = pd.read_csv('/Data_2/Jun/Adipocytes/references/representative-isoforms.txt', sep='\t', names=['Gene_ID','Transcript_ID','a','b','c']).drop(columns=['a','b','c'])

DB = pd.merge(repre, P_gtf, on=['Gene_ID','Transcript_ID'])
DB = DB.rename(columns = {'Gene_ID':'gene_id',
                          'Gene_name':'gene_name',
                          'Transcript_ID':'transcript_id'})
DB = DB[['transcript_id',"5UTR_length","CDS_length","3UTR_length","gene_id","gene_name"]]
DB = DB.set_index('transcript_id').T.to_dict('list')

prefix = 'adi_'
#data file
folder = '/Data_2/Jun/Adipocytes/tr-aln/novaseq/'
days = ['D0','D4','D8']
reps = ['a','b','c']
filetype = ".rep.bam"

#%%
#constants
Ribosome_shift_factor = [-21,-10] #CDS range

result = {'Frame':["F1", "F2", "F3"]}
for day in days:
    for rep in reps:
        SP = day+rep
        result[SP] = [0]*3
        direc = folder + SP + filetype
        f = pysam.AlignmentFile(direc,'rb')
        for line in f:
            line = line.tostring(f)
            Aligned = line.split("\t")
            
            T_id = Aligned[2]
            if DB.get(T_id)==None: continue
            P = DB[T_id]
            
            start_pos = int(Aligned[3]) - int(P[0]+1)
            stop_pos = int(Aligned[3]) - int(P[0]+P[1]-2)
            
            if start_pos >= Ribosome_shift_factor[0] and stop_pos <= Ribosome_shift_factor[1] : #CDS part
                result[SP][ (start_pos - Ribosome_shift_factor[0]) % 3 ] += 1
        print(SP+' done')
            
result = pd.DataFrame(result)
result.to_csv(prefix+'3ntP.tsv', index=False, sep='\t')    