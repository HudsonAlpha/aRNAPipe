# -*- coding: utf-8 -*-
import sys
import os

# LIBRARY USED TO SUBMIT JOBS: 'LSF' (IBM LSF WORKLOAD MANAGER), 'LOCAL' (SEQUENTIAL RUN ON SINGLE MACHINE) OR 'OTHER' FOR LIBRARY ADAPTED TO OTHER WORKLOAD MANAGERS
mode = "LSF"
# PATH TO THE FOLDER "genomes_processed" WHERE THE DIFFERENT GENOME BUILDS ARE STORED
path_db = "/gpfs/gpfs1/myerslab/reference/genomes/rnaseq_pipeline/bedpipe/"
# PATHS TO BINARIES USED BY aRNApipe
path_fastqc   = "/gpfs/gpfs1/software/fastqc/fastqc"
path_kallisto = "/gpfs/gpfs1/software/kallisto-0.42.4/kallisto"
path_star     = "/gpfs/gpfs1/software/STAR_2.4.2a/bin/Linux_x86_64_static/STAR"
path_htseq    = "/gpfs/gpfs1/software/HTSeq-0.5.3/bin/htseq-count"
path_picard   = "/gpfs/gpfs1/software/picard-tools-1.88"
path_samtools = "/gpfs/gpfs1/software/samtools-1.2/bin/samtools"
path_varscan  = "/gpfs/gpfs1/software/varscan/VarScan.v2.3.6.jar"
path_gtf2gp   = "/gpfs/gpfs1/myerslab/reference/genomes/rnaseq_pipeline/bin/gtfToGenePred" # last version
path_gatk     = "/gpfs/gpfs1/software/GATK-3.5/GenomeAnalysisTK.jar"
path_cutadapt = "/gpfs/gpfs1/software/python2.7/bin/cutadapt"
path_trimgalore = "/gpfs/gpfs1/myerslab/reference/genomes/rnaseq_pipeline/bin/trim_galore"
path_starfusion = "/gpfs/gpfs1/software/STAR_2.4.2a/STAR-Fusion-0.1.1/STAR-Fusion"

# STAR options: The keys of this dict are used in the project config files to use the referenced STAR arguments
star_options = {"default": "",
                "encode":  "--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000",
                "fusion":  "--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 10 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 200000 --alignMatesGapMax 200000"}

# ENVIRONMENT VARIABLES TO CHANGE (ONLY FOR THE PIPELINE EXECUTION)
environment = {"JAVA_HOME":      ["/usr/java/jdk1.8.0_60/","add"],
               "PYTHONPATH":     ["/gpfs/gpfs1/software/HTSeq-0.5.3/lib/python","overwrite"],
               "PATH":           ["/gpfs/gpfs1/software/Python-2.7.2/bin","add"],
               "LD_LIBRARY_PATH":["/gpfs/gpfs1/software/gcc-4.8.2/usr/lib64","add"],
               "PERL5LIB"       :["/gpfs/gpfs1/software/perl-modules/lib/perl5/5.10.1:/gpfs/gpfs1/software/perl-modules/lib/perl5/5.10.1/lib64/perl5","add"]}

# ANNOTATIONS AND FULL PATH TO THE PIPELINE BASE DIRECTORY
path_genome = path_db + "/genomes_processed/#LABEL/STAR_genome"
path_index = path_db + "/genomes_processed/#LABEL/kallisto.idx"
path_annotation = path_db + "/genomes_processed/#LABEL/genesets.gtf"
path_fasta = path_db + "/genomes_processed/#LABEL/genome.fa"
annots = [path_db + "/genomes_processed/#LABEL/genesets.refFlat",
          path_db + "/genomes_processed/#LABEL/refFlats/protein_coding.refFlat",
          path_db + "/genomes_processed/#LABEL/refFlats/rRNA.refFlat"]
nannots = ["general","protein_coding","ribosomal"]

# FILES REQUIRED BY GATK
gatk_multithread = {"RTC": 4, "BR": 4, "PR": 4}
annots_gatk = {"g1k_v37": ["dbsnp_138.b37.vcf", ["1000G_phase1.indels.b37.vcf", "Mills_and_1000G_gold_standard.indels.b37.vcf"]]}
