# -*- coding: utf-8 -*-
import sys
import os

# LIBRARY USED TO SUBMIT JOBS:
# - 'LSF' FOR IBM LSF WORKLOAD MANAGER (it uses 'sys_LSF.py')
# - 'LOCAL' FOR SEQUENTIAL RUN ON SINGLE MACHINE (it uses 'sys_single.py')
# - 'OTHER' FOR LIBRARIES ADAPTED TO OTHER WORKLOAD MANAGERS (it uses 'sys_OTHER.py')
mode = "LSF"

# PATH TO THE FOLDER "genomes_processed" WHERE THE DIFFERENT GENOME BUILDS ARE STORED
path_db = "/gpfs/gpfs1/myerslab/reference/genomes/rnaseq_pipeline/"

# FULL PATHS TO BINARIES USED BY aRNApipe (users must change these values to match
# the current locations of the binaries used by aRNApipe in their system).
path_fastqc   = "/gpfs/gpfs1/software/fastqc/fastqc" #aut
path_kallisto = "/gpfs/gpfs1/software/kallisto-0.42.4/kallisto" #aut
path_star     = "/gpfs/gpfs1/software/STAR_2.5.2b/bin/Linux_x86_64_static/STAR"
path_htseq    = "/gpfs/gpfs1/software/HTSeq-0.5.3/bin/htseq-count" # no
path_picard   = "/gpfs/gpfs1/software/picard-tools-1.88" #aut
path_samtools = "/gpfs/gpfs1/software/samtools-1.2/bin/samtools" #aut
path_varscan  = "/gpfs/gpfs1/software/varscan/VarScan.v2.3.6.jar" #aut
path_gtf2gp   = "/gpfs/gpfs1/myerslab/reference/genomes/rnaseq_pipeline/bin/gtfToGenePred" #aut
path_gatk     = "/gpfs/gpfs1/software/GATK-3.5/GenomeAnalysisTK.jar" # no
path_cutadapt = "/gpfs/gpfs1/software/python2.7/bin/cutadapt" # no
path_trimgalore = "/gpfs/gpfs1/myerslab/reference/genomes/rnaseq_pipeline/bin/trim_galore" # aut
path_starfusion = "/gpfs/gpfs1/software/STAR_2.5.2b/STAR-Fusion/STAR-Fusion"
path_jsplice = "/gpfs/gpfs1/software/jSplice-1.0.1/" #aut

# STAR options (users can add their own options):
# The keys of this dict are used in the project config files to use the
# referenced STAR arguments within the corresponding dictionary values
star_options = {"default": "",
                "encode":  "--outFilterType BySJout --outFilterMultimapNmax 20 --alignSJoverhangMin 8 --alignSJDBoverhangMin 1 --outFilterMismatchNmax 999 --outFilterMismatchNoverLmax 0.04 --alignIntronMin 20 --alignIntronMax 1000000 --alignMatesGapMax 1000000"}

# ENVIRONMENT VARIABLES:
# The following system environment variables are changed to add or overwrite
# their current values.
environment = {"JAVA_HOME":      ["/usr/java/jdk1.8.0_60/","add"],
               "PYTHONPATH":     ["/gpfs/gpfs1/software/HTSeq-0.5.3/lib/python","overwrite"],
               "PATH":           ["/gpfs/gpfs1/software/Python-2.7.2/bin","add"],
               "PATH":           ["/gpfs/gpfs1/software/bedtools2-2.20.0/bin","add"],
               "PATH":           ["/gpfs/gpfs1/software/samtools-1.2/bin","add"],
               "LD_LIBRARY_PATH":["/gpfs/gpfs1/software/gcc-4.8.2/usr/lib64","add"],
               "PERL5LIB"       :["/gpfs/gpfs1/software/perl-modules/lib/perl5/5.10.1:/gpfs/gpfs1/software/perl-modules/lib/perl5/5.10.1/lib64/perl5","add"]}

# ANNOTATIONS AND FULL PATH TO THE PIPELINE BASE DIRECTORY (do not change)
path_genome = path_db + "/genomes_processed/#LABEL/STAR_genome"
path_index = path_db + "/genomes_processed/#LABEL/kallisto.idx"
path_annotation = path_db + "/genomes_processed/#LABEL/genesets.gtf"
path_star_fusion = path_db + "/genomes_processed/#LABEL/star-fusion"
path_fasta = path_db + "/genomes_processed/#LABEL/genome.fa"
annots = [path_db + "/genomes_processed/#LABEL/genesets.refFlat",
          path_db + "/genomes_processed/#LABEL/refFlats/protein_coding.refFlat",
          path_db + "/genomes_processed/#LABEL/refFlats/rRNA.refFlat"]
nannots = ["general","protein_coding","ribosomal"]

# ARGUMENTS REQUIRED BY GATK
# Check GATK documentation for RTC, BR and PR values:
# http://gatkforums.broadinstitute.org/dsde/discussion/1975/how-can-i-use-parallelism-to-make-gatk-tools-run-faster
gatk_multithread = {"RTC": 4, "BR": 4, "PR": 4}
# key: Genome version label
# values: (list)
# - First element: SNPs
# - Second element: List of known indels files used by BaseRecalibrator, RealignerTargetCreator and IndelRealigner
# https://software.broadinstitute.org/gatk/guide/article?id=1247
annots_gatk = {}
annots_gatk["g1k_v37"] = ["dbsnp_138.b37.vcf", ["1000G_phase1.indels.b37.vcf", "Mills_and_1000G_gold_standard.indels.b37.vcf"]]
