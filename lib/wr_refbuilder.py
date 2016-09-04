import optparse
import time
import os
import config as config
import refbuild as refbuild

##########################################################
## Parsing command line arguments
##########################################################
desc = "aRNApipe: Reference builder"
parser = optparse.OptionParser(description = desc)
parser.add_option("-L", "--label", dest = "label",default = "", help = "[Required] Identifier label for the genome.")
parser.add_option("-p", "--path",  dest = "path", default = "", help = "[Required] Absolute path to the directory where the folder 'genomes_processed' is located.")
parser.add_option("-f", "--fasta", dest = "fasta",default = "", help = "Path to the uncompressed genome fasta file ('.fa' or '.fasta').")
parser.add_option("-c", "--cdna",  dest = "cdna", default = "", help = "Path to the cDNA fasta file (accepts '.gz').")
parser.add_option("-g", "--gtf",   dest = "gtf",  default = "", help = "Path to the GTF gene set file ('.gtf').")
parser.add_option("-n", "--ncpu",  dest = "n",    default = "8",help = "Number of threads that STAR will use to generate the reference genome (default=8).")

(opt, args) = parser.parse_args()

(opt, args) = parser.parse_args()

root_path = opt.path
file_fasta = opt.fasta
file_cdna = opt.cdna
file_gtf = opt.gtf
genome_label = opt.label
nprocs_star  = int(opt.n)

print "Options: "
print "- Genome label:     " + genome_label
print "- Pipeline path:    " + root_path
print "- Input GTF file:   " + file_gtf
print "- Input FASTA file: " + file_fasta
print "- Input cDNA file:  " + file_cdna

## CHECKS RAW AND PROCESSED ANNOTATIONS FOLDER
print "> Checking processed genome directories..."
if not os.path.exists(root_path + "/genomes_processed"):
    os.mkdir(root_path + "/genomes_processed")
    print "  Directory for processed genomes created."
else:
    print "  Directory for processed genomes already exists."

## CHECKS GENOME LABEL
print "> Checking if genome label is already installed..."
if os.path.exists(root_path + "/genomes_processed/installed_genomes.txt"):
    f = open(root_path + "/genomes_processed/installed_genomes.txt", 'r')
    for i in f:
        if i.strip("\n") == genome_label:
            exit("A genome with the same label has already been installed: " + genome_label)
    f.close()

## CHECKS INPUT FILES
k = 0
print "> Checking input fasta/cdna/gtf files..."
if file_fasta != "":
    if not os.path.exists(file_fasta):
        exit("Some of the input GTF/CDNA/FASTA files not found.")
else:
    print " FASTA file not provided."
    k += 1
if file_cdna != "":
    if not os.path.exists(file_cdna):
        exit("Some of the input GTF/CDNA/FASTA files not found.")
else:
    print " cDNA file not provided."
    k += 1
if file_gtf != "":
    if not os.path.exists(file_gtf):
        exit("Some of the input GTF/CDNA/FASTA files not found.")
else:
    print " GTF file not provided."
    k += 1

## CREATES OUTPUT DIRECTORY
print "> Creating output directory for the current genome version in the processed genomes folder..."
if os.path.exists(root_path + "/genomes_processed/" + genome_label):
    os.system("rm -r " + root_path + "/genomes_processed/" + genome_label)
os.mkdir(root_path + "/genomes_processed/" + genome_label)
if not os.path.exists(root_path + "/genomes_processed/" + genome_label + "/log"):
    os.mkdir(root_path + "/genomes_processed/" + genome_label + "/log")
if not os.path.exists(root_path + "/genomes_processed/" + genome_label + "/temp"):
    os.mkdir(root_path + "/genomes_processed/" + genome_label + "/temp")
path_logs = root_path + "/genomes_processed/" + genome_label + "/log/"
path_temp = root_path + "/genomes_processed/" + genome_label

## COPYING ORIGINAL FILES
if file_fasta != "":
    print "> Making a copy of the FASTA file to the processed genome directory: 'genome.fa'"
    os.system("cp " + file_fasta + " " + root_path + "/genomes_processed/" + genome_label + "/genome.fa")
if file_cdna != "":
    print "> Making a copy of the cDNA file to the processed genome directory: 'transcriptome.fa'"
    os.system("cp " + file_cdna + " " + root_path + "/genomes_processed/" + genome_label + "/transcriptome.fa")
if file_gtf != "":
    print "> Making a copy of the GTF file to the processed genome directory: 'genesets.gtf'"
    os.system("cp " + file_gtf + " " + root_path + "/genomes_processed/" + genome_label + "/genesets.gtf")

## GENERATE ANNOTATION
if file_gtf != "":
    print "> Generating exon, transcript and gene annotation from GTF..."
    res = refbuild.annotate_gtf(root_path + "/genomes_processed/" + genome_label + "/genesets.gtf")
    if res == 0:
        exit("Error generating annotations...")

## GENERATE INDEXES FOR GATK
if file_fasta != "":
    print "> Generating index and dict files..."
    res = refbuild.fasta2dict(root_path + "/genomes_processed/" + genome_label + "/genome.fa", config.path_samtools, config.path_picard)
    if res == 0:
        exit("Error generating dict and index files...")

## STAR: CREATE GENOME FOR STAR ANALYSIS
if file_fasta != "":
    print "> Creation of reference genome for STAR analysis..."
    print "  - Path: " + root_path + "/genomes_processed/" + genome_label + "/STAR_genome"
    if not os.path.exists(root_path + "/genomes_processed/" + genome_label + "/STAR_genome"):
        os.mkdir(root_path + "/genomes_processed/" + genome_label + "/STAR_genome")
    if file_gtf != "":
        print "  - Including GTF annotation"
        command  = config.path_star + " --runThreadN " + str(nprocs_star) + " --runMode genomeGenerate  --genomeDir #GD --genomeFastaFiles #FASTA --sjdbGTFfile #GTF &>/dev/null"
    else:
        print "  - GTF file not included/provided"
        command  = config.path_star + " --runThreadN " + str(nprocs_star) + " --runMode genomeGenerate  --genomeDir #GD --genomeFastaFiles #FASTA &>/dev/null"
    command  = command.replace("#GD", root_path + "/genomes_processed/" + genome_label + "/STAR_genome")
    command  = command.replace("#FASTA", file_fasta).replace("#GTF", file_gtf)
    os.system(command)

## KALLISTO: CREATE TRANSCRIPTOME INDEX FOR KALLISTO ANALYSIS
if file_cdna != "":
    print "> Creation of transcriptome index for Kallisto analysis: 'kallisto.idx'"
    command  = config.path_kallisto + "  index -i #OUTPUT #INPUT &>/dev/null"
    command  = command.replace("#OUTPUT", root_path + "/genomes_processed/" + genome_label + "/kallisto.idx")
    command  = command.replace("#INPUT", file_cdna)
    os.system(command)

## REFFLAT FORMAT
if file_gtf != "":
    print "> Converting GTF file to refFlat format..."
    print "> Splitting GTF file by feature type and converting to refFlat format..."
    res = refbuild.refflat(config.path_gtf2gp, root_path + "/genomes_processed/" + genome_label)
    if res == 0:
        exit("Error generating reFflat files...")

## annotate
print "> Adding genome key to installed genomes..."
out = open(root_path + "/genomes_processed/installed_genomes.txt",'a')
print >> out, genome_label + "\t" + time.strftime("%y%m%d_%H%M%S") + "\t" + opt.fasta.split("/")[-1] + "\t" + opt.cdna.split("/")[-1] + "\t" + opt.gtf.split("/")[-1]
out.close()
