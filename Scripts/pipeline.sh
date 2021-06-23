#!/bin/bash
#BSUB -o out.FQ2NA_Run_010621.%J
#BSUB -e err.FQ2NA_Run_010621.%J
#BSUB -W 48:00
#BSUB -n 2
#BSUB -M 7000
#BSUB -R "span[ptile=1]"
#BSUB -J "FQ2NA_Run_010621"



### Files with sample names present in each patient's folder
# $pat_data/Samplelist.txt;
# $pat_data/Samplelist_tumour.txt;
# $pat_data/Samplelist_normal.txt;

#RNA
data_RNA="/gpfs/projects/vhio76/vhio76640/Intership_Laura/RNA"
quantiseq_path="/gpfs/projects/vhio76/vhio76640/Intership_Laura/RNA"

while read pat; do

	pat_data="$data_RNA/$pat.data"
	rm -f $pat_data/RNA_seq_files_${pat}.txt

	while read sample; do

		fastp -i $pat_data/${pat}-${sample}_RNA_1.fq -I $pat_data/${pat}-${sample}_RNA_2.fq -o $pat_data/${pat}-${sample}_RNA_1.fq.gz -O $pat_data/${pat}-${sample}_RNA_2.fq.gz

		printf "${pat}-${sample}\t$pat_data/${pat}-${sample}_RNA_1.fq.gz\t$pat_data/${pat}-${sample}_RNA_2.fq.gz\n" >> $pat_data/RNA_seq_files_${pat}.txt

	done < $pat_data/Samplelist.txt; #Contains sample names for that patient

  printf "$pat prepatation done!\t"

	$quantiseq_path/quanTIseq_pipeline.sh --inputfile=$pat_data/RNA_seq_files_${pat}.txt --outputdir=$pat_data --tumor=TRUE --prefix=Sample --rawcounts=TRUE --threads=2

	printf "$pat quanTIseq DONE!\n"

done < $data/Patientlist.txt; #Contains all the patients names

printf "all patients done\n"

echo "Quantification of the Tumor Immune contexture";


# WES



date

module load gcc
module load fastp
module load bwa
module load samtools

data="/gpfs/projects/vhio76/vhio76640/Intership_Laura/WES/Sample_323" ## change this to the data directory

reference="/gpfs/projects/vhio76/vhio76640/Intership_Laura/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fna"

samtools faidx /gpfs/projects/vhio76/vhio76640/Intership_Laura/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fna

bwa index $reference

while read pat; do

	pat_data="$data/$pat.data"

	while read sample; do

		fastp -i $pat_data/${sample}_1.fq -I $pat_data/${sample}_2.fq -o $pat_data/${sample}_1.fq.gz -O $pat_data/${sample}_2.fq.gz

		bwa mem -aM -R "@RG\tID:${sample}\tPU:${sample}.1\tSM:${sample}\tPL:ILLUMINA\tPI:150" $reference $pat_data/${sample}_1.fq.gz $pat_data/${sample}_2.fq.gz | samtools view -bS - > $pat_data/${sample}.bam


		samtools index $pat_data/${sample}.bam
		samtools sort $pat_data/${sample}.bam > $pat_data/${sample}_S.bam
		samtools index $pat_data/${sample}_S.bam

       ### UNCOMMENT TO REMOVE THE SORTED FILE AFTER USING IT
       # rm $pat_data/${sample}.bam

	done < $pat_data/Samplelist.txt; #Contains sample names for that patient

    printf "$pat\t"

 done < $data/Patientlist.txt;

printf "all patients done\n"

echo "Mapping FASTQs DONE!";


#Optitype

module purge
module load intel/14.0.2 impi/4.1.3.049 gcc/5.1.0 OPENSSL/1.1.1c PERL/5.20.2 PYTHON/3.5.1-INTEL NETMHCPAN/4.0 RAZERS3 SAMTOOLS/1.2 GLPK/4.65 HDF5/1.8.10

while read pat; do

	pat_data="$data/$pat.data"

	while read sample; do

		python /apps/OPTIPYTE/1.3.2/OptiTypePipeline.py -i $pat_data/${sample}_1.fq.gz $pat_data/${sample}_2.fq.gz --dna -c config.ini --verbose -o $pat_data/.

	done < $pat_data/Samplelist.txt;

  printf "$pat\t"

done < $data/Patientlist.txt;

printf "all patients done\n"

echo "OPTITYPE DONE!";

 
 module purge
 
 module load gcc
 module load intel
 module load spark
 module load gatk
 module load samtools

 #SAMTOOLS GATK
 
 while read pat; do
 
 	pat_data="$data/$pat.data"
 
 	while read sample; do
 
 		samtools rmdup $pat_data/${sample}_S.bam $pat_data/${sample}_RD.bam
 
 		gatk MarkDuplicatesSpark --input $pat_data/${sample}_RD.bam --output $pat_data/${sample}_001.bam --tmp-dir $pat_data/tmp
 
         ### UNCOMMENT TO REMOVE THE SORTED FILE AFTER USING IT
         # rm $pat_data/${sample}_S.bam
         # rm $pat_data/${sample}_RD.bam
 
 	done < $pat_data/Samplelist.txt;
 
     printf "$pat\t"
 
 done < $data/Patientlist.txt;
 
 printf "all patients done\n"
 
 echo "REMOVING DUPLICATES DONE!";
 
 module purge
 
 module load java/1.8.0u102
 module load gatk
 
 reference="/gpfs/projects/vhio76/vhio76640/Intership_Laura/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fna"
 ##Double check this file!
 knownsites="/gpfs/projects/vhio76/vhio76640/NeoAnt/Data_PL/common_all_20180418.vcf.gz"
 
 gatk CreateSequenceDictionary -R $reference
 
 while read pat; do
 	pat_data="$data/$pat.data"
 	while read sample; do
 		gatk BaseRecalibrator -I Sample_323.data/323-007_001.bam -R $reference --known-sites $knownsites -O Sample_323.data/323-007_recal_data.table
 		gatk ApplyBQSR -R $reference -I Sample_323.data/323-007_001.bam  --bqsr-recal-file Sample_323.data/323-007_recal_data.table -O Sample_323.data/323-007_002.bam 
 
 	done < $pat_data/Samplelist.txt;
    printf "$pat\t"
 
 done < $data/Patientlist.txt;
 
 printf "all patients done\n"
 echo "RECALIBRATION DONE!";
 
############
############
############
# UNTIL HERE ALL THE SAMPLES ARE PROCESSED THE SAME WAY
#      FROM NOW ON, THE *NORMAL* SAMPLE IS USED DIFFERENTLY FROM THE OTHERS
############
############
############

odule purge

odule load spark
odule load gatk


reference="/gpfs/projects/vhio76/vhio76640/Intership_Laura/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fna"
germline_gnomad="/gpfs/projects/vhio76/vhio76640/NeoAnt/Data_PL/somatic-hg38_af-only-gnomad.hg38.vcf";
while read pat; do

	pat_data="$data/$pat.data"
	normal_sample=$(head -1 $pat_data/Samplelist_normal.txt);

	while read sample; do

		gatk Mutect2 -R $reference -I $pat_data/${sample}_002.bam -I $pat_data/${normal_sample}_002.bam -normal ${normal_sample} -tumor ${sample} --germline-resource $germline_gnomad --f1r2-tar-gz $pat_data/${sample}_f1r2.tar.gz -O $pat_data/${sample}_somatic.vcf.gz

		gatk LearnReadOrientationModel -I $pat_data/${sample}_f1r2.tar.gz -O $pat_data/${sample}_read-orientation-model.tar.gz

	done < $pat_data/Samplelist_tumour.txt; 

    printf "$pat\t"

done < $data/Patientlist.txt;

printf "all patients done\n"
echo "Mutect2 DONE!";

somatic_exac="/gpfs/projects/vhio76/vhio76640/NeoAnt/Data_PL/somatic-hg38_small_exac_common_3.hg38.vcf.gz";

while read pat; do
	pat_data="$data/$pat.data"

	while read sample; do

		# Tabulates pileup metrics for inferring contamination
		gatk GetPileupSummaries -I $pat_data/${sample}_002.bam -V $somatic_exac -L $somatic_exac -O $pat_data/${sample}_pileups.table

	done < $pat_data/Samplelist.txt;
    printf "$pat\t"

done < $data/Patientlist.txt;

printf "all patients done\n"
echo "PILEUP CONTAMINATION DONE!";


reference="/gpfs/projects/vhio76/vhio76640/Intership_Laura/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fna"
germline_gnomad="/gpfs/projects/vhio76/vhio76640/NeoAnt/Data_PL/somatic-hg38_af-only-gnomad.hg38.vcf";
while read pat; do

	pat_data="$data/$pat.data"
	normal_sample=$(head -1 $pat_data/Samplelist_normal.txt);

	while read sample; do

		# Calculate the fraction of reads coming from cross-sample contamination
		gatk CalculateContamination -I $pat_data/${sample}_pileups.table -matched $pat_data/${normal_sample}_pileups.table --tumor-segmentation $pat_data/${sample}_segments.table -O $pat_data/${sample}_contamination.table

		# Filter somatic SNVs and indels called by Mutect2
		gatk FilterMutectCalls -V $pat_data/${sample}_somatic.vcf.gz -R $reference --tumor-segmentation $pat_data/${sample}_segments.table --contamination-table $pat_data/${sample}_contamination.table --ob-priors $pat_data/${sample}_read-orientation-model.tar.gz -O $pat_data/${sample}_somatic_filtered.vcf

	done < $pat_data/Samplelist_tumour.txt; #Contains tumour sample names for that patient

    printf "$pat\t"

done < $data/Patientlist.txt;

printf "all patients done\n"

echo "COMPUTE CONTAMINATION & FILTERING MUTECT DONE!";







## Select variants that pass the filters and Funcotator
reference="/gpfs/projects/vhio76/vhio76640/Intership_Laura/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fna"
## Double check path!
data_sources_path="/gpfs/scratch/vhio76/vhio76640/NeoAnt_Pipeline/Funcotator/dataSourcesFolder/"

while read pat; do
	pat_data="$data/$pat.data"

	while read sample; do

		grep -E "#|PASS" $pat_data/${sample}_somatic_filtered.vcf > $pat_data/${sample}_somatic_filtered_PASS.vcf

		#Funcotator (FUNCtional annOTATOR) analyzes given variants for their function (as retrieved from a set of data sources) and produces the analysis in a specified output file.
		gatk Funcotator -R $reference -V $pat_data/${sample}_somatic_filtered_PASS.vcf -O $pat_data/${sample}_somatic_annotated.vcf --output-file-format VCF --data-sources-path $data_sources_path --ref-version hg38

	done < $pat_data/Samplelist_tumour.txt; 

    printf "$pat\t"

done < $data/Patientlist.txt;

printf "all patients done\n"

echo "SELECT VARIANTS & FUNCOTATOR DONE!";


#Sequenza

module purge
module load gcc/8.4.0 openmpi/4.0.2 MKL/2017.4 R/4.0.4
module load intel/2017.4 impi/2017.4 MKL/2017.4 gcc/8.4.0 OPENSSL/1.1.1c PYTHON/3.7.4_pip


## Add sequenza
data="/gpfs/projects/vhio76/vhio76640/Intership_Laura/WES/Sample_323"; ### remember to change!
date


reference="/gpfs/projects/vhio76/vhio76640/Intership_Laura/GRCh38/GCA_000001405.15_GRCh38_full_analysis_set.fna.gz"
wig_file="/gpfs/scratch/vhio76/vhio76640/NeoAnt_Pipeline/Sequenza/hg38.gc20Base.wig.gz"; # add path to wig file !!

path_script="$data"; # path to R script


while read pat; do
	pat_data="$data/$pat.data"
	normal_sample=$(head -1 $pat_data/Samplelist_normal.txt);

	input_normal_bam="$pat_data/${normal_sample}_S.bam"

	while read sample; do

		input_tumour_bam="$pat_data/${sample}_S.bam"

		# bam2seqz creating a file per chr
		while read chr; do
			date
			echo sequenza-utils bam2seqz -n $input_normal_bam -t $input_tumour_bam --fasta $reference -gc $wig_file -o $pat_data/${sample}_${chr}_DELETE.seqz.gz -C $chr
			date
		done < $data/chr_list.txt;
        printf "Chromosomes files CREATED!\t"



		# bam2seqz merging the chr files
		sequenza-utils seqz_merge -o $pat_data/${sample}_a_DELETE.seqz.gz -1 $pat_data/${sample}_chr1_DELETE.seqz.gz -2 $pat_data/${sample}_chr2_DELETE.seqz.gz
		while read chr; do
			date
			sequenza-utils seqz_merge -o $pat_data/${sample}_b_DELETE.seqz.gz -1 $pat_data/${sample}_a_DELETE.seqz.gz -2 $pat_data/${sample}_${chr}_DELETE.seqz.gz
			mv $pat_data/${sample}_b_DELETE.seqz.gz $pat_data/${sample}_a_DELETE.seqz.gz
			date
		done < $data/chr_list_minus1_2.txt;

		mv $pat_data/${sample}_a_DELETE.seqz.gz $pat_data/${sample}.seqz.gz

		# remove unwanted intermediate files
		rm $pat_data/${sample}*_DELETE.seqz.gz
        printf "Chromosomes MERGED!\t"


		# sequenza STEP 2
		sequenza-utils seqz_binning --seqz $pat_data/${sample}.seqz.gz -w 20 -o $pat_data/${sample}_small.seqz.gz
        printf "Binning DONE!\t"

		# sequenza STEP 3

		Rscript $path_script/sequenza_results.R $pat_data/${sample}_small.seqz.gz ${sample}
        printf "Rscript DONE!\n"

	done < $pat_data/Samplelist_tumour.txt;

    echo "$pat DONE!\n"

done < $data/Patientlist.txt;

printf "all patients done\n"

echo "SEQUENZA DONE!";




##Pyclone

module purge

module load gcc/latest openssl/1.0.1l python/2.7.14 llvm/8.0.1

path_script="$data"; # path to Python script

while read pat; do

	pat_data="$data/$pat.data"
	mkdir -p $pat_data/tmp;

    	export NUMBA_CACHE_DIR=$pat_data/tmp;

			# generate copy number variants files
			## ** REQUIRES pat-sample_segments.txt files from Sequenza **
			python $path_script/run_pyclone.py $pat_data $pat_data/Samplelist_tumour.txt $pat;

			# run PyClone with the output genereated by the previous script
			PyClone run_analysis_pipeline --in_files $pat_data/*_CNV.tsv --working_dir $pat_data --density pyclone_beta_binomial --num_iters 10000 --prior total_copy_number


    echo "$pat DONE!"

done < $data/Patientlist.txt;

printf "all patients done\n"

echo "PyClone DONE!";




# NeopredPipe

while read pat; do
	pat_data="$data/$pat.data"

	while read sample; do

		NeoPredPipe -I $pat_data -H $pat_data/${sample}.HLA.tsv -o $pat_data -n ${sample} -m

	done < $pat_data/Samplelist.txt;

  echo "$pat DONE!\n"

done < $data/Patientlist.txt;

printf "all patients done\n"

echo "NeoPredPipe DONE!";


