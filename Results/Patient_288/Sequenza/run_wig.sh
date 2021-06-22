#!/bin/bash
#BSUB -o out.Sequenza_wig_20_160621.%J
#BSUB -e err.Sequenza_wig_20_160621.%J
#BSUB -W 48:00
#BSUB -n 2
#BSUB -M 7000
#BSUB -R "span[ptile=1]"
#BSUB -J "Sequenza_wig_20_160621"

module purge && module load intel/2017.4 impi/2017.4 MKL/2017.4 gcc/8.4.0 OPENSSL/1.1.1c PYTHON/3.7.4_pip

date 

sequenza-utils gc_wiggle -w 20 --fasta GCA_000001405.15_GRCh38_full_analysis_set.fna -o GRCh38.gc20Base.wig.gz

date
