# -*- coding: utf-8 -*-

import os
import shutil
import sys
import time
import optparse
import lib.config as config
import lib.vcrparser as vcrparser
# Dynamics load of the workload manager library depending on 'config.mode' value set in lib/config.py
if config.mode == "LSF":
    import lib.sys_LSF as manager
elif config.mode == "LOCAL":
    import lib.sys_single as manager
else:
    import lib.sys_OTHER as manager

##########################################################
## Parsing command line arguments
##########################################################
desc = "aRNApipe: RNA-seq framework"
parser = optparse.OptionParser(description = desc)
parser.add_option("-m", "--mode", dest = "m", default = "",         help = "[Required] Write mode: 'new'/'update'/'kill'/'progress'/'skeleton'/'genomes'. 'new' all the libraries will be processed and previously processed results, if exist, will be removed; 'update' only those libraries that have not been previously processed will be analyzed; 'kill' stops a current execution of the given project by killing all its processes; 'progress' shows the progress of a running execution; 'skeleton' creates a skeleton for a project using the absolute path to the project folder given by -p; 'genomes' displays the genome versions currently available.")
parser.add_option("-p", "--project_folder", dest = "folder", default = "", help = "[Required] Absolute path to the project folder. This folder will be used to store the results and must contain the files 'samples.list' and 'config.txt'.")
parser.add_option("-w", "--walltime", dest = "wt", default = "200:00",     help = "Optional: Wall time for the main job (defaults to 200:00).")
parser.add_option("-q", "--queue", dest = "q", default = "priority",       help = "Optional: Queue that will be used by the main job (defaults to 'priority').")
(opt, args) = parser.parse_args()
path_aRNApipe = os.path.dirname(sys.argv[0])

##########################################################
## Genomes mode: Shows available genome builds and exits
##########################################################
if opt.m == "genomes":
    db = path_aRNApipe+ "/../genomes_processed/installed_genomes.txt"
    if os.path.exists(db):
        f = open(db, 'r')
        print "Available genome builds:"
        for i in f:
            i = i.strip("\n").split("\t")
            if len(i) >= 2:
                print "- Key: " + i[0] + " (installed " + i[1] +"): " + "/".join(i[2:])
        f.close()
        exit("All Done!")
    else:
        exit("Database file not found.")

#############################################################
## Skeleton mode: Builds skeleton for a new project and exits
#############################################################
if opt.m == "skeleton":
    if not os.path.exists(opt.folder):
        try:
            os.mkdir(opt.folder)
            shutil.copy(path_aRNApipe + "/template/config.txt", opt.folder + "/config.txt")
            shutil.copy(path_aRNApipe + "/template/samples.list", opt.folder + "/samples.list")
            print "Project skeleton created."
        except:
            exit("Unspecified error during skeleton generation.")
    else:
        print "Project folder already exists"
    exit("All Done!")

##########################################################
## Checking required arguments (path and mode)
##########################################################
if (opt.folder == "") or (opt.m == ""):
    exit("Parameters not provided. Check --help.")
if not opt.folder.startswith("/"):
    exit("An absolute path must be provided (started with '/').")
if opt.folder.endswith("/"):
    opt.folder = opt.folder[0:-1]
opt.folder    = opt.folder.split("/")
opt.path_base = "/".join(opt.folder[0:-1])
opt.folder    = opt.folder[-1]
complete_path = opt.path_base + "/" + opt.folder

##########################################################
## Check for user confirmation in update, kill or new modes
##########################################################
if opt.m == "kill":
    t = raw_input("Kill the current analysis (y/n)?")
    if t != "y":
        exit("--> Exit")
elif opt.m == "update":
    print "When running 'update' be sure that no previuos processes on the same project are running."
    t = raw_input("Continue (y/n)?")
    if t != "y":
        exit("--> Exit")
elif opt.m == "new":
    print "Mode 'new'. All the results data in the project folder, if exists, will be removed."
    t = raw_input("Continue (y/n)?")
    if t != "y":
        exit("--> Exit")

#######################################################################
## Removing previous results, temp folder and log files if mode = 'new'
#######################################################################
if opt.m == "new":
    n = os.listdir(complete_path)
    for i in n:
        if i.startswith("results_"):
            nom = complete_path + "/" + i
            os.system("rm -r " + nom.replace("//","/"))
        elif i.startswith("aRNApipe_"):
            nom = complete_path + "/" + i
            os.system("rm " + nom.replace("//","/"))
        elif i=="temp":
            nom = complete_path + "/" + i
            os.system("rm -r " + nom.replace("//","/"))
        elif i=="pid.txt":
            nom = complete_path + "/" + i
            os.system("rm " + nom.replace("//","/"))
        elif i=="logs":
            nom = complete_path + "/" + i
            os.system("rm -r " + nom.replace("//","/"))
elif opt.m == "update":
    n = os.listdir(complete_path)
    for i in n:
        if i.startswith("results_"):
            if os.path.exists(complete_path + "/" + i + "/samples_ko.txt"):
                os.remove(complete_path + "/" + i + "/samples_ko.txt")

##########################################################
## Submiting the main process if mode 'new' or 'update'
##########################################################
print "aRNApipe:"
print "- Input: " + complete_path
print "- Mode:  " + opt.m
## Kill current analysis if mode 'kill'
if   opt.m == "kill":
    n = vcrparser.project_kill(opt.path_base, opt.folder)
## Check progress if mode 'progress'
elif opt.m == "progress":
    n = vcrparser.project_process(opt.path_base, opt.folder)
## Submiting the main process if mode 'new' or 'update'
elif opt.m in ["update", "new"]:
    timestamp = time.strftime("%y%m%d_%H%M%S")
    if not os.path.exists(complete_path):
        os.mkdir(complete_path)
    if not os.path.exists(complete_path + "/logs"):
        os.mkdir(complete_path + "/logs")
    ## submittint job
    vcr_args  = " -f " + opt.folder + " -b " + opt.path_base + " -m " + opt.m + " -t " + timestamp
    if opt.m == "update":
        out = open(complete_path + "/logs/aRNApipe.log", 'a')
        print >> out, "###########################################"
        print >> out, "## Update: " + timestamp
    else:
        out = open(complete_path + "/logs/aRNApipe.log", 'w')
        print >> out, "###########################################"
        print >> out, "## New: " + timestamp
    print >> out, "###########################################"
    out.close()
    if not os.path.exists(complete_path + "/temp"):
        os.mkdir(complete_path + "/temp")
    ce      = vcrparser.change_environment(config.environment)
    uds = manager.submit_job(opt.wt, "1", opt.q, complete_path + "/logs/aRNApipe_cluster_" + timestamp + ".log", opt.folder, "python " + path_aRNApipe + "/lib/wr_aRNApipe.py" + vcr_args + " >> " + complete_path + "/logs/aRNApipe.log", 0, complete_path, "")
    print "Main process submitted (jid=" + uds + ")"
    out = open(complete_path + "/pid.txt",'w')
    print >> out, "aRNApipe\t" + uds + "\t" + timestamp
    out.close()
else:
    exit("Selected mode not valid.")
