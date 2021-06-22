## STEP 1 & 2

import io
import sys
import os
import pandas as pd
import numpy as np
import re


# READ COMMAND LINE ARGUMENTS
path_to_vcfs = sys.argv[1] # provide the full path -> WITHOUT THE LAST "/"
sample_list_file = sys.argv[2] # provide the full path to the file
patient_id = sys.argv[3] # only the patient identifier
normal_sample = sys.argv[4] # full sample name pat-sam e.g.(288-024)

print(path_to_vcfs)
print(sample_list_file)
print(patient_id)
print(normal_sample)


## Create list of the different samples for that patient
with open(sample_list_file, 'r') as f:
    different_samples = [l.strip() for l in f]
    print(different_samples)

# at this point if all the input has been read properly the program should work.



def read_vcf(path, name, name_normal):
    with open(path, 'r') as f:
        lines = [l for l in f if not l.startswith('#')]
    return pd.read_csv(
        io.StringIO(''.join(lines)),
        dtype={'#CHROM': str, 'POS': str, 'ID': str, 'REF': str, 'ALT': str,
               'QUAL': str, 'FILTER': str, 'INFO': str, name: str, name_normal: str},
        sep='\t'
    ).rename(columns={'#CHROM': 'CHROM'})



def clean_vcf(vcf,clean_vcf, name, name_normal):
    VCF = read_vcf(vcf, name, name_normal)
    TUMOUR = (VCF[name])
    NORMAL = (VCF[name_normal])
    INFO = (VCF['INFO'])
    CHROM = (VCF['CHROM'])
    POS = (VCF['POS'])

    listGENE = []
    listGN_CH_POS_ID = []
    listREF_COUNTS = []
    listVAR_COUNTS = []

    #print(INFO)

    for index, row in INFO.items():
        #print(row)
        #t = row.split('|')
        t = re.split('\||;|\[|\]',row)
        #print(t[4])

        listGENE.append(t[4]) # we select gene name
        listGN_CH_POS_ID.append(t[4] + '_' + CHROM[index] + '_' + POS[index]) # we build the mutation id
        # break

    #print("info done")

    for index, row in TUMOUR.items():
        t = row.split(':')
        s = t[1].split(',')
        listREF_COUNTS.append(s[0])

    for index, row in TUMOUR.items():
        t = row.split(':')
        s = t[1].split(',')
        listVAR_COUNTS.append(s[1])

    d = ('mutation_id', 'ref_counts', 'var_counts', 'normal_cn', 'minor_cn', 'major_cn')
    VCF_Filtered = pd.DataFrame(columns=d)

    VCF_Filtered['mutation_id'] = pd.Series(listGN_CH_POS_ID)
    VCF_Filtered['ref_counts'] = pd.Series(listREF_COUNTS)
    VCF_Filtered['var_counts'] = pd.Series(listVAR_COUNTS)
    VCF_Filtered['normal_cn'] = '2'
    VCF_Filtered['minor_cn'] = '1'
    VCF_Filtered['major_cn'] = '1'

    VCF_Filtered.to_csv(clean_vcf, sep='\t', index=False)
    print("VCF cleaned and saved")


def vcf2list(name, name_normal):
    #input
    vcf_filename = path_to_vcfs + "/" + name + "_somatic_annotated.vcf"

    #output
    filename = path_to_vcfs + "/" + name + ".tsv"

    clean_vcf(vcf_filename, filename, name, name_normal)
    DF = pd.read_csv(filename, sep="\t")
    set0 = set(DF["mutation_id"].to_numpy())

    return set0



# Read all the files and concatenate the mutations
setofmutations = set()

for sample in different_samples:
    name_of_file = patient_id + '_' + sample

    setofmutations = setofmutations.union(vcf2list(name_of_file, normal_sample))

list_all = list(setofmutations)
# print(list_all[:20])

DF_all = pd.DataFrame(list_all)
DF_all.to_csv(path_to_vcfs + "/" + patient_id + "_AllMut.tsv", index=False, header= False, sep="\t")

# print(DF_all.columns)
print("DF_all saved to file")


#It could be useful to remove the mutations without HUGO symbol
#Edit the output to get three columns "chr" "start" & "end" (start + 1) this will be the bed file for next step
all_mut_file = path_to_vcfs + "/" + patient_id + "_AllMut"

#query = "sed 's/_/\t/1' " + all_mut_file + ".tsv | rev | sed 's/_/\t/1' | rev | awk 'BEGIN{FS=OFS=\"\t\"}{print $2\"\t\"$3\"\t\"$3+1\"\t\"$1}' | sed 's/ /\t/g' > " + all_mut_file + ".txt"
#os.system(query)

query2 = "sed 's/_/\t/1' 288_AllMut.tsv | rev | sed 's/_/\t/1' | rev | awk 'BEGIN{FS=OFS=\"\t\"}{print $2\"\t\"$3\"\t\"$3+1}' | paste - " + all_mut_file + ".tsv > " + all_mut_file + ".txt"
os.system(query2)

#new_df = DF_all[["chr", "start"]]
#new_df["end"] = new_df["start"] + 1
#new_df.to_csv(path_to_vcfs + "/" + patient_id + "_AllMut.txt", index=False, sep="\t")
print("DF_all txt saved to file")



## STEP 3
# REMOVE UNWANTED MUTATION TYPES
os.system("sed -i '/random/d' *.tsv")
os.system("sed -i '/chrY/d' *.tsv")
os.system("sed -i '/chrUn/d' *.tsv")
os.system("sed -i '/chrM/d' *.tsv")
os.system("sed -i '/alt/d' *.tsv")
os.system("sed -i '/\./d' *.tsv")

print("messy chromosomes removed")



# Clean A
#From AllMut file create regions_AllMut.bed


## STEP 4
regionsfile_patient = path_to_vcfs + "/" + patient_id + "_AllMut.txt"

for sample in different_samples:
    name_of_file = patient_id + '_' + sample


    # we first sort the bam file
    bam_file_input = path_to_vcfs + "/" + name_of_file + "_002.bam"
    bam_file_out_sorted = path_to_vcfs + "/" + name_of_file + "_002_S.bam" # name a bit different from before

    command_sort = "samtools sort " + bam_file_input + " > " + bam_file_out_sorted
    os.system(command_sort)
    print("SAMTOOLS done!")


    counts_file = path_to_vcfs + "/" + name_of_file + "_Counts.tsv"

    command = "bedtools coverage -counts -a " + regionsfile_patient + " -b " + bam_file_out_sorted + " > " + counts_file
    os.system(command)
    print("BEDTOOLS done!")

    ## UNCOMMENT TO REMOVE THE SORTED FILE AFTER USING IT
    # command_remove = "rm " + bam_file_out_sorted
    # os.system(command_remove)




## STEP 5
def Merge_TSV(Sample_TSV, Sample_Counts, path_OUT):
    TSV_Sam = pd.read_csv(Sample_TSV, sep="\t")
    print(TSV_Sam.columns)
    Cou_Sam = pd.read_csv(Sample_Counts, sep="\t")
    print(Cou_Sam.columns)
    New_TSV_File = pd.DataFrame()
    New_TSV_File.columns = ["mutation_id", "total"]
    New_TSV_File = TSV_Sam.merge(Cou_Sam, left_on=['mutation_id'], right_on=['mutation_id'], how='right')
    New_TSV_File.to_csv(path_OUT, index=False, sep="\t")


for sample in different_samples:
    name_of_file = patient_id + '_' + sample

    # inputs
    tsv_vcf = path_to_vcfs + "/" + name_of_file + ".tsv"
    counts_file = path_to_vcfs + "/" + name_of_file + "_Counts.tsv"

    # output
    out_all_mut_file = path_to_vcfs + "/" + name_of_file + "_AllMut.tsv"

    Merge_TSV(tsv_vcf, counts_file, out_all_mut_file)





## STEP 6

def Add_CNV_2_TSV(path_CN, path_TSV, path_OUT):
    CN_File = pd.read_csv(path_CN, sep="\t")
    TSV_File = pd.read_csv(path_TSV, sep="\t")

    listChr = []
    listPos = []

    INFO = (TSV_File['mutation_id'])

    for index, row in INFO.items():
        t = row.split('_')
        listChr.append(t[1])
        listPos.append(t[2])


    TSV_File['Chromosome'] = pd.Series(listChr)
    TSV_File['Pos'] = pd.Series(listPos)

    chDic = {}
    posCol = []
    posCol2 = []
    posCol3 = []

    for idx,ID2 in enumerate(CN_File['chromosome']):
        if ID2 not in chDic:
            chDic[ID2] = []
        chDic[ID2].append([CN_File['start.pos'][idx], CN_File['end.pos'][idx], CN_File['CNt'][idx], CN_File['B'][idx], CN_File['A'][idx]])

    for idx,ID2 in enumerate(TSV_File['Chromosome']):
        if ID2 in chDic:
            match = 2
            pos = TSV_File['Pos'][idx]
            for idx2, ID3 in enumerate(chDic[ID2]):
                try:
                    if int(pos) >= int(ID3[0]) and int(pos) <= int(ID3[1]):
                        match = 2
                        match2 = ID3[3]
                        match3 = ID3[4]
                except:
                    match = 2
                    match2 = 1
                    match3 = 1
            posCol.append(match)
            posCol2.append(match2)
            posCol3.append(match3)
        else:
            posCol.append(2)
            posCol2.append(1)
            posCol3.append(1)

    TSV_File['normal_cn'] = pd.Series(posCol)
    TSV_File['minor_cn'] = pd.Series(posCol2)
    TSV_File['major_cn'] = pd.Series(posCol3)
    TSV_File = TSV_File.drop(columns=['Chromosome', 'Pos'])
    listHugo = []
    Filter = (TSV_File["mutation_id"])

    for index, row in Filter.items():
        t = row.split('_')
        listHugo.append(t[0])

    TSV_File['Hugo_Symb'] = pd.Series(listHugo)
    TSV_File = TSV_File.replace(r'^\s*$', np.nan, regex=True)
    TSV_File = TSV_File.dropna(subset=["Hugo_Symb"])
    TSV_File = TSV_File.drop(columns=['Hugo_Symb'])
    TSV_File.to_csv(path_OUT, index=False, sep="\t")



for sample in different_samples:
    name_of_file = patient_id + '_' + sample

    # inputs
    segments_file = path_to_vcfs + "/" + name_of_file + "_segments.txt" # Sequenza output
    all_mut_file = path_to_vcfs + "/" + name_of_file + "_AllMut.tsv"

    # output
    out_CNV = path_to_vcfs + "/" + name_of_file + "_CNV.tsv"

    Add_CNV_2_TSV(segments_file, all_mut_file, out_CNV)



    
