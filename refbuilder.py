import sys
import optparse
import os

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
parser.add_option("-w", "--wt",    dest = "wt",   default = "200:00", help = "Wall time (default=200:00).")
parser.add_option("-s", "--span",  dest = "span",   default = "no", help = "Span 1 host (yes/no, default no).")
(opt, args) = parser.parse_args()

# Span one host?
if opt.span == "no":
    g = ""
else:
    g = "-R span[hosts=1] "

# Build arguments for 'wr_refbuilder.py'
print "> Creating output directory for the current genome version in the processed genomes folder..."
if os.path.exists(opt.path + "/genomes_processed/" + opt.label):
    os.system("rm -r " + opt.path + "/genomes_processed/" + opt.label)
os.mkdir(opt.path + "/genomes_processed/" + opt.label)
if not os.path.exists(opt.path + "/genomes_processed/" + opt.label + "/log"):
    os.mkdir(opt.path + "/genomes_processed/" + opt.label + "/log")
if not os.path.exists(opt.path + "/genomes_processed/" + opt.label + "/temp"):
    os.mkdir(opt.path + "/genomes_processed/" + opt.label + "/temp")

vargs = "-L " + opt.label + " -p " + opt.path + " -f " + opt.fasta + " -c " + opt.cdna + " -g " + opt.gtf + " -n " + opt.n
bsub_1 = "bsub " + g + "-q normal -J " + opt.label + " -n " + opt.n +" -W " + opt.wt + " -o " + opt.path + "/genomes_processed/" + opt.label + '/' + opt.label + "_cluster.log"
bsub_2 = " 'python " + os.path.dirname(sys.argv[0]) + "/lib/wr_refbuilder.py " + vargs + " > " + opt.path + "/genomes_processed/" + opt.label + '/' + opt.label + ".log"
os.system(bsub_1 + bsub_2)
