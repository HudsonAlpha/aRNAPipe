# -*- coding: utf-8 -*-
import os
import time
import optparse
import config
import vcrparser
import programs
import shutil
if config.mode == "LSF":
    import sys_LSF as manager
elif config.mode == "LOCAL":
    import sys_single as manager
else:
    import sys_OTHER as manager

########################################################################
## HELP OPTIONS AND DOCUMENTATION
desc = "aRNApipe: RNA-seq framework"
parser = optparse.OptionParser(description=desc)
parser.add_option("-f", "--project_folder", dest = "folder",    default = "",  help = "")
parser.add_option("-m", "--write",          dest = "m",         default = "0", help = "")
parser.add_option("-b", "--path_base",      dest = "path_base", default = "",  help = "")
parser.add_option("-t", "--timestamp",      dest = "timestamp", default = "",  help = "")
########################################################################
(opt, args) = parser.parse_args()
timestamp = opt.timestamp
print "###########################################"
print "### aRNApipe: RNA-seq framework    #########"
print "###########################################"
print "> Analysis started: " + timestamp
print "> Parsing options..."
## PATH BASE AND PROJECT FOLDER CHECK
[path_base, folder] = vcrparser.check_args(opt.path_base, opt.folder)
## PARSE CONFIGURATION FILE
opt.samples = path_base + folder + "/samples.list"
opt.config  = path_base + folder + "/config.txt"
[opt.config, var]   = vcrparser.config_file(opt.config, path_base, folder, config)
shutil.copy(opt.samples, opt.samples.replace("samples.list","logs/" + timestamp + "_samples.list"))
shutil.copy(opt.config, opt.config.replace("config.txt","logs/" + timestamp + "_config.txt"))
## CHECK GENOME BUILD
config.path_genome = config.path_genome.replace("#LABEL", var["genome_build"])
if not os.path.exists(config.path_genome):
    exit("path_genome not found. Genome build " + var["genome_build"] + " missing or incomplete.")
config.path_index = config.path_index.replace("#LABEL", var["genome_build"])
if not os.path.exists(config.path_index):
    exit("path_index not found. Genome build " + var["genome_build"] + " missing or incomplete.")
config.path_annotation = config.path_annotation.replace("#LABEL", var["genome_build"])
if not os.path.exists(config.path_annotation):
    exit("path_annotation not found. Genome build " + var["genome_build"] + " missing or incomplete.")
for i in range(3):
    config.annots[i] = config.annots[i].replace("#LABEL", var["genome_build"])
    if not os.path.exists(config.annots[i]):
        exit("annots not found. Genome build " + var["genome_build"] + " missing or incomplete.")
## VERBOSE OPTIONS
print "> GLOBAL VARIABLES:"
print "  - Project folder:  " + folder
print "  - Base path:       " + path_base
print "  - Config file:     " + opt.config
print "  - Genome build:    " + var["genome_build"]
print "  - Run mode:        " + opt.m
print "  - Stranded RNA:    " + var["strandedness"]
print "> CLUSTER VARIABLES:"
print "  - Queue:           " + var["q"]
print "  - Walltime:        " + var["wt"]
print "> ANALYSIS:"
print "  - TrimGalore:      " + var["trimgalore"]
print "  - FastQC:          " + var["fastqc"]
print "  - STAR:            " + var["star"]
print "  - STAR-Fusion:     " + var["star-fusion"]
print "  - Picard QC:       " + var["picard"]
print "  - HTseq (gene):    " + var["htseq-gene"]
print "  - HTseq (exon):    " + var["htseq-exon"]
print "  - Kallisto:        " + var["kallisto"]
print "  - Varscan:         " + var["varscan"]
print "> INDIVIDUAL ANALYSIS SETTINGS:"
print "  - TrimGalore args: " + var["trimgal_args"]
print "  - STAR arguments:  " + var["star_args"]
print "  - STAR 2-pass:     " + var["star2pass"]
print "  - STAR output-fus: " + var["starfusion"]
print "  - Kall bootstraps: " + var["kalboot"]
print "  - VARSCAN args:    " + var["varscan_args"]
print "  - GATK args:       " + var["gatk_args"]
print "  - HTseqGene mode:  " + var["htseq-gene-mode"]
print "  - HTseqExon mode:  " + var["htseq-exon-mode"]


samples = vcrparser.get_samples(path_base, folder, opt.samples)
##########################################################
## Creates 'temp' folder for temporary managing files
##########################################################
n = os.listdir(opt.path_base + "/" + opt.folder)
if not ("temp" in n):
    os.mkdir(path_base + "/" + folder + "/temp")
else:
    out = open(path_base + "/" + folder + "/temp/pids.txt", 'w')
    out.close()
    out = open(path_base + "/" + folder + "/temp/pids_scheduler.txt", 'w')
    out.close()
##########################################################
## Strand-specific protocol and HTseq mode
##########################################################
if var["strandedness"] == "yes":
    var["strandedness"] = " --stranded=yes"
elif var["strandedness"] == "reverse":
    var["strandedness"] = " --stranded=reverse"
elif var["strandedness"] == "no":
    var["strandedness"] = " --stranded=no"
else:
    exit("Error: Strandedness value not correct (yes/no/reverse).")
if (not (var["htseq-gene-mode"] in ["union", "intersection-strict", "intersection-nonempty"])) or (not (var["htseq-exon-mode"] in ["union", "intersection-strict", "intersection-nonempty"])):
    exit("Error: HTseq mode value not correct (union/intersection-strict/intersection-nonempty).")

##########################################################
## Check STAR options
##########################################################
if int(var["star"].split("/")[0]) > 0:
    if var["star_args"] == "own":
            var["star_args"] = var["star_args_own"]
    else:
        if config.star_options.has_key(var["star_args"]):
            var["star_args"] = config.star_options[var["star_args"]]
        else:
            exit("star_args key not found in 'config.py' file.")
    if var["starfusion"] == "own":
        var["starfusion"] = var["starfusion_own"]
    else:
        var["starfusion"] = "--chimSegmentMin 12 --chimJunctionOverhangMin 12"
    if var["star2pass"] == "yes":
        var["star2pass"] = "--twopassMode Basic"
    else:
        var["star2pass"] = ""
    args = dict()
    for i in ["star2pass", "star_args", "starfusion"]:
        if var[i] != "":
            j = var[i].split(" ")
            if (len(j) % 2) != 0:
                exit(i + " not properly formatted.")
            for k in range(0, len(j), 2):
                if not j[k].startswith("--"):
                    exit(i + " not properly formatted.")
                args[j[k]] = j[k + 1]
    star_params = ""
    if len(args) > 0:
        for i, j in args.iteritems():
            star_params = star_params + " " + i + " " + j
##########################################################
## Other checks
##########################################################
if int(var["gatk"].split("/")[0]) > 0:
    if not config.annots_gatk.has_key(var["genome_build"]):
        exit("GATK annotation files are not available for this genome build: " + var["genome_build"])
if int(var["star-fusion"].split("/")[0]) > 0:
    if len(samples[samples.keys()[0]]) == 2:
        exit("Star-Fusion requires paired-end reads.")

##########################################################
## Starts analysis
##########################################################
procs = list()
tg = False
## Pre-alignment analysis
if int(var["trimgalore"].split("/")[0]) > 0:
    tg = True
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "trimgalore", opt.m)
    if len(samples_v) > 0:
        uds_tg, logs_tg = programs.trimgalore(timestamp, path_base, folder, samples_v, var["trimgalore"], var["wt"], var["q"], var["trimgal_args"])
        w = vcrparser.job_wait(logs_tg, 20)
        procs.append(logs_tg)
if int(var["fastqc"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "fastqc", opt.m)
    if len(samples_v) > 0:
        uds_qc, logs_qc = programs.fastqc(timestamp, path_base, folder, samples_v, var["fastqc"], var["wt"], var["q"], tg)
        procs.append(logs_qc)
if int(var["kallisto"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "kallisto", opt.m)
    if len(samples_v) > 0:
        uds_kal, logs_kal = programs.kallisto(timestamp, path_base, folder, samples_v, config.path_index, var["kalboot"], var["kallisto"], var["wt"], var["q"], tg)
        procs.append(logs_kal)
if int(var["star"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "star", opt.m)
    if len(samples_v) > 0:
        uds_star, logs_star = programs.star(timestamp, path_base, folder, samples_v, var["star"], var["wt"], var["q"], config.path_genome, star_params, tg)
        w = vcrparser.job_wait(logs_star, 20)
        procs.append(logs_star)
if int(var["star-fusion"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "star-fusion", opt.m)
    if len(samples_v) > 0:
        uds_sf, logs_sf  = programs.starfusion(timestamp, path_base, folder, samples_v, var["star-fusion"], var["wt"], var["q"], var["genome_build"])
        procs.append(logs_sf)
if int(var["picard"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "picard", opt.m)
    if len(samples_v) > 0:
        uds_qca, logs_qca  = programs.picardqc(timestamp, path_base, folder, samples_v, var["picard"], var["wt"], var["q"], config.annots, var["strandedness"])
        procs.append(logs_qca)
if int(var["htseq-gene"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "htseq-gene", opt.m)
    if len(samples_v) > 0:
        uds_htseqG, logs_htseqG = programs.htseq(timestamp, path_base, folder, samples_v, config.path_annotation, var["htseq-gene"], var["wt"], var["q"], "gene", var["strandedness"], var["htseq-gene-mode"])
        procs.append(logs_htseqG)
if int(var["htseq-exon"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "htseq-exon", opt.m)
    if len(samples_v) > 0:
        uds_htseqE, logs_htseqE = programs.htseq(timestamp, path_base, folder, samples_v, config.path_annotation, var["htseq-exon"], var["wt"], var["q"], "exon", var["strandedness"],var["htseq-exon-mode"])
        procs.append(logs_htseqE)
if (int(var["varscan"].split("/")[0]) > 0) or (int(var["gatk"].split("/")[0]) > 0):
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "sam2sortbam", opt.m)
    if len(samples_v) > 0:
        uds_sam2sortbam, logs_sam2sortbam = programs.sam2sortbam(timestamp, path_base, folder, samples_v, var["varscan"], var["wt"], var["q"])
        procs.append(logs_sam2sortbam)
        w = vcrparser.job_wait(logs_sam2sortbam, 20)
if int(var["gatk"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "gatk", opt.m)
    if len(samples_v) > 0:
        uds_gatk, logs_gatk = programs.gatk(timestamp, path_base, folder, samples_v, var["gatk"], var["wt"], var["q"], var["genome_build"], var["gatk_args"])
        procs.append(logs_gatk)
if int(var["varscan"].split("/")[0]) > 0:
    samples_v, stats = vcrparser.check_samples(samples, path_base, folder, "varscan", opt.m)
    if len(samples_v) > 0:
        uds_varscan, logs_varscan = programs.varscan(timestamp, path_base, folder, samples_v, var["varscan"], var["wt"], var["q"], var["genome_build"], var["varscan_args"])
        procs.append(logs_varscan)

if len(procs) > 0:
    for proc in procs:
        w = vcrparser.job_wait(proc, 10)

timestamp = time.strftime("%y%m%d_%H%M%S")
print "> Analysis finished: " + timestamp

