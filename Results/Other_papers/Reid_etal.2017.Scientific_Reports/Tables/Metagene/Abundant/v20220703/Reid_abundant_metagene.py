import pandas as pd
import pysam

col_list = ["Transcript_ID","Gene_ID","Gene_name","5UTR_length","CDS_length","3UTR_length"]
P_gtf = pd.read_csv('/Data_1/Daehwa/python/GTF_Parsing/v2.2/gencode.vM27.annotation.gtf/v20220414/Processed_gtf.tsv', sep='\t', usecols=col_list)
repre = pd.read_csv('/Data_2/Daehwa/Adipocyte/Other_papers/Reid_etal.2017.Scientific_Reports/Alignment/references/representative-isoforms.txt', sep='\t', names=['Gene_ID','Transcript_ID','a','b','c']).drop(columns=['a','b','c'])
abundant = pd.read_csv('/Data_2/Daehwa/Adipocyte/Other_papers/Reid_etal.2017.Scientific_Reports/Analysis/Gene_lists/Abundant_genes/v20220703/Reid_RPF_top100_abundant_genelist.txt', names=['gene_id'])

DB = pd.merge(repre, P_gtf, on=['Gene_ID','Transcript_ID'])
DB = DB.rename(columns = {'Gene_ID':'gene_id',
                          'Gene_name':'gene_name',
                          'Transcript_ID':'transcript_id'})
DB = pd.merge(abundant, DB, on='gene_id')
DB = DB[['transcript_id',"5UTR_length","CDS_length","3UTR_length","gene_id","gene_name"]]
DB = DB.set_index('transcript_id').T.to_dict('list')

prefix = 'Reid_'
#data file
folder = '/Data_2/Daehwa/Adipocyte/Other_papers/Reid_etal.2017.Scientific_Reports/Alignment/tr-aln/RPF/'
days = ['BD0','BD5','WD0','WD5','BAT','WAT']
reps = ['a','b']
filetype = ".rep.bam"

#%%
for day in days:
    for rep in reps:
        SP = prefix+day+rep
        result_start = {}
        result_stop = {}
        result_start['transcript_id'] = [str(i) for i in range(-100,201)]
        result_stop['transcript_id'] = [str(i) for i in range(-200,101)]
        
        for T_id in DB:
            result_start[T_id] = [0]*len(range(-100,201))
            result_stop[T_id]  = [0]*len(range(-200,101))
        
        direc = folder + SP + filetype
        f = pysam.AlignmentFile(direc,'rb')
        for line in f:
            line = line.tostring(f)
            Aligned = line.split("\t")
            
            T_id = Aligned[2]
            
            if DB.get(T_id)==None: continue
            P = DB[T_id]

            start_pos = int( int(Aligned[3])-(P[0]+1) )
            stop_pos = int( int(Aligned[3])-(P[0]+P[1]-2) )
            if -100 <= start_pos <= 200 :
                result_start[T_id][start_pos+100] += 1
            if -200 <= stop_pos <= 100 :
                result_stop[T_id][stop_pos+200] += 1
        
        print (SP+' done')
        
        result_start = pd.DataFrame(result_start).T
        result_start.to_csv(SP+'_abndt_metagene-start.tsv', header=False, sep='\t')
        result_stop = pd.DataFrame(result_stop).T
        result_stop.to_csv(SP+'_abndt_metagene-stop.tsv', header=False, sep='\t')
        