# -*- coding: utf-8 -*-
import os
import config
import time

if config.mode == "LSF":
    import sys_LSF as manager
elif config.mode == "LOCAL":
    import sys_single as manager
else:
    import sys_OTHER as manager

def change_environment(env):
    if len(env) > 0:
        for var, val in env.iteritems():
            if val[1] == "overwrite":
                os.environ[var] = val[0]
            else:
                if os.environ.has_key(var):
                    os.environ[var] = val[0] + ":" + os.environ[var]
                else:
                    os.environ[var] = val[0]

def project_process(path_base, folder):
    samples = get_samples(path_base, folder, path_base + "/" + folder + "/samples.list")
    # Check main process
    print "## MAIN PROCESS ###########################"
    try:
        f = open(path_base + "/" + folder + "/pid.txt", 'r')
        i = f.readline().strip("\n").split("\t")[1]
        f.close()
        k = manager.job_status(i)
        if k == 1:
            st = "DONE"
        elif k == 0:
            st = "RUN"
        else:
            st = "ERROR"
        print "- Main process (" + i + ") status: " + st
    except:
        print "- Main process not found or already finished"
    # Check subprocesses
    print "## SUBPROCESSES ###########################"
    pids = dict()
    try:
        f = open(path_base + "/" + folder + "/temp/pids.txt", 'r')
        for i in f:
            i = i.strip("\n").split("\t")
            pids[i[0]] = [i[1].split("|"),i[2].split("|")]
        f.close()
    except:
        print "- No subprocesses file found"
    f = open(path_base + "/" + folder + "/config.txt", 'r')
    config = dict()
    for i in f:
        if not i.startswith("%"):
            i = i.strip("\n").split("\t")
            if i[0] in ["trimgalore", "fastqc", "kallisto", "star", "star-fusion", "picard", "htseq-gene", "htseq-exon", "varscan", "gatk"]:
                i[1] = i[1].split("/")[0]
                if i[1] != "0":
                    config[i[0]] = i[1]
    if config.has_key("varscan") or config.has_key("gatk"):
        config["sam2sortbam"] = 1
    if len(config) > 0:
        for pg in ["trimgalore", "fastqc", "kallisto", "star", "star-fusion", "picard", "htseq-gene", "htseq-exon", "sam2sortbam", "varscan", "gatk"]:
            if config.has_key(pg):
                print "Process:  " + pg
                if not pids.has_key(pg):
                    print "- Already done or waiting for previous module output"
                else:
                    pid = pids[pg]
                    print "- ID:       " + "|".join(pid[0])
                    n = list()
                    for i in pid[1]:
                        k = manager.job_status(i)
                        if k == 1:
                            n.append("DONE")
                        elif k == 0:
                            n.append("RUN")
                        else:
                            n.append("ERROR")
                    print "- Status:   " +  "|".join(n)
                    samples_v, stats = check_samples(samples, path_base, folder, pg, "update")
                    sok = str(round(100 * float(stats[1])/float(stats[0]),2))
                    sko = str(round(100 * float(stats[2])/float(stats[0]),2))
                    pending = str(round(100 * float(stats[0]-stats[1]-stats[2])/float(stats[0]),2))
                    print "- Progress: " + sok + "% succeeded / " + sko + "% exited / " + pending + "% pending"


#############################################################################
# MAIN FUNCTION FOR KILLING ALL THE PROCESSES RELATED TO A PROJECT RUN
#############################################################################
def project_kill(path_base, folder):
    print "Main process:"
    uds = ""
    try:
        f = open(path_base + "/" + folder + "/pid.txt",'r')
        uds = f.readline().strip("\n").split("\t")[1]
        f.close()
        print "- Killing main process ("+uds+")"
        manager.job_kill(uds)
    except:
        print "- No main process to kill. Already finished (" + uds + ")"
    print "Submodule processes:"
    try:
        f = open(path_base + "/" + folder + "/temp/pids.txt", 'r')
        for i in f:
            i = i.strip("\n").split("\t")
            if len(i) > 1:
                j = i[1].split("|")
                for jj in j:
                    print "- Killing process " + i[0] + " ("+jj+")"
                    manager.job_kill(jj)
        f.close()
    except:
        print "No submodule processes to kill"
    return 1


#############################################################################
## PARSES DE CONFIGURATION FILE
#############################################################################
def config_file(config, path_base, folder, paths):
    # config: path to configuration file
    # path_base: Absolute path to the location where the project folder has been created
    # folder: Name of the project folder located in 'path_base'
    f = open(config, 'r')
    var = dict()
    for i in f:
       if not i.startswith("%"):
           i = i.rstrip().split("\t")
           if len(i) > 1:
               var[i[0]] = i[1]
    f.close()
    return [config, var]


def job_wait(path, secs):
## WAITS UNTIL A JOB ENDS SUCCESSFULLY (IF ERROR -> ABORTS) BASED ON THE LOG FILES CREATED BY THE RESOURCE MANAGER
## Input:
##   path: String with the path to the log files corresponding to the jobs that will be checked (i.e. single job: '/some_path/job1.log', job array: '/some_path/job1.log|/some_path/job2.log')
##   secs: Number of seconds between subsequent checks
## Returns:
##   1: Once all the jobs have finished
## If an error is detected -> exit()
    if "|" in path:
        jobarray = path.split("|")
    else:
        jobarray = [path]
    for path in jobarray:
        while 1:
            rt = manager.job_status(path)
            if rt == 1:
                break
            elif rt == -1:
                exit("Error on: " + path)
            else:
                time.sleep(secs)
    return 1

#############################################################################
## PARSES DE SAMPLES FILE
#############################################################################
def get_samples_nocheck(path_base, folder, samplefile):
    # SCAN FOR FASTQ SAMPLE FILES
    try:
        f = open(samplefile, 'r')
    except:
        exit("Error: Samples file not found")
    samples = dict()
    # CHECK COLUMN HEADERS AND SINGLE-END/PAIRED-END DATA
    i = f.readline().strip("\n").split("\t")
    idx = [-1, -1, -1]
    for j in range(len(i)):
        if i[j] == "SampleID":
            idx[0] = j
        elif i[j] == "FASTQ_1":
            idx[1] = j
        elif i[j] == "FASTQ_2":
            idx[2] = j
        elif i[j] == "FASTQ":
            idx[1] = j
    # 'SampleID' AND 'FASTQ' COLUMNS ARE REQUIRED
    if (idx[0] < 0) or (idx[1] < 0):
        exit("Error: Samples file headers must contain 'SampleID' and 'FASTQ' columns for single-end or 'SampleID', 'FASTQ_1' and 'FASTQ_2' for paired-end data")
    # PARSE SAMPLE DATA
    errors = dict({"ID duplication errors":[],"Missing input files":[], "Empty input files":[]})
    for i in f:
        i = i.strip("\n").split("\t")
        if len(i) > 1:
            # PAIRED-END (TWO FASTQ FILES PROVIDED PER SAMPLE
            if idx[2] >= 0:
                # CHECKS FILE EXTENSION OK (*.fastq or *.fastq.gz)
                if (i[idx[1]].endswith("fastq") and i[idx[2]].endswith("fastq")) or (i[idx[1]].endswith("fastq.gz") and i[idx[2]].endswith("fastq.gz")):
                    # NO DUPLICATE IDS
                    if samples.has_key(i[idx[0]]):
                        errors["ID duplication errors"].append(i[idx[0]])
                    else:
                        for ifile in range(1,3):
                            samples[i[idx[0]]] = [i[idx[1]], i[idx[2]], 0, 0]
                else:
                    exit("Error: Input sample files must be '.fastq' or '.fastq.gz'")
            # SINGLE-END (ONE FASTQ FILES PROVIDED PER SAMPLE
            else:
                # CHECKS FILE EXTENSION OK (*.fastq or *.fastq.gz)
                if i[idx[1]].endswith("fastq") or i[idx[1]].endswith("fastq.gz"):
                    # NO DUPLICATE IDS
                    if samples.has_key(i[idx[0]]):
                        errors["ID duplication errors"].append(i[idx[0]])
                    else:
                        for ifile in range(1,2):
                            samples[i[idx[0]]] = [i[idx[1]], 0]
                else:
                    exit("Error: Input sample files must be '.fastq' or '.fastq.gz'")
    if len(samples) == 0:
        exit("Error: No available samples identified.")
    r = 0
    for i,j in errors.iteritems():
        if len(j) > 0:
            r += 1
            print i + ":"
            for k in j:
                print "- " + k
    if r > 0:
        exit("Samples file errors detected")
    return samples


def get_samples(path_base, folder, samplefile):
    # SCAN FOR FASTQ SAMPLE FILES
    try:
        f = open(samplefile, 'r')
    except:
        exit("Error: Samples file not found")
    samples = dict()
    # CHECK COLUMN HEADERS AND SINGLE-END/PAIRED-END DATA
    i = f.readline().strip("\n").split("\t")
    idx = [-1, -1, -1]
    for j in range(len(i)):
        if i[j] == "SampleID":
            idx[0] = j
        elif i[j] == "FASTQ_1":
            idx[1] = j
        elif i[j] == "FASTQ_2":
            idx[2] = j
        elif i[j] == "FASTQ":
            idx[1] = j
    # 'SampleID' AND 'FASTQ' COLUMNS ARE REQUIRED
    if (idx[0] < 0) or (idx[1] < 0):
        exit("Error: Samples file headers must contain 'SampleID' and 'FASTQ' columns for single-end or 'SampleID', 'FASTQ_1' and 'FASTQ_2' for paired-end data")
    # PARSE SAMPLE DATA
    errors = dict({"ID duplication errors":[],"Missing input files":[], "Empty input files":[]})
    for i in f:
        i = i.strip("\n").split("\t")
        if len(i) > 1:
            # PAIRED-END (TWO FASTQ FILES PROVIDED PER SAMPLE
            if idx[2] >= 0:
                # CHECKS FILE EXTENSION OK (*.fastq or *.fastq.gz)
                if (i[idx[1]].endswith("fastq") and i[idx[2]].endswith("fastq")) or (i[idx[1]].endswith("fastq.gz") and i[idx[2]].endswith("fastq.gz")):
                    # NO DUPLICATE IDS
                    if samples.has_key(i[idx[0]]):
                        errors["ID duplication errors"].append(i[idx[0]])
                    else:
                        for ifile in range(1,3):
                            # CHECKS IF INPUT FILES EXIST
                            try:
                                print i[idx[ifile]]
                                f = open(i[idx[ifile]], 'r')
                                f.close()
                                S = os.stat(i[idx[ifile]]).st_size
                                if S > 0:
                                    samples[i[idx[0]]] = [i[idx[1]], i[idx[2]], os.stat(i[idx[1]]).st_size, os.stat(i[idx[2]]).st_size]
                                else:
                                    errors["Empty input files"].append(i[idx[ifile]])
                            except:
                                errors["Missing input files"].append(i[idx[ifile]])
                else:
                    exit("Error: Input sample files must be '.fastq' or '.fastq.gz'")
            # SINGLE-END (ONE FASTQ FILES PROVIDED PER SAMPLE
            else:
                # CHECKS FILE EXTENSION OK (*.fastq or *.fastq.gz)
                if i[idx[1]].endswith("fastq") or i[idx[1]].endswith("fastq.gz"):
                    # NO DUPLICATE IDS
                    if samples.has_key(i[idx[0]]):
                        errors["ID duplication errors"].append(i[idx[0]])
                    else:
                        for ifile in range(1,2):
                            # CHECKS IF INPUT FILES EXIST
                            try:
                                f = open(i[idx[ifile]],'r')
                                f.close()
                                S = os.stat(i[idx[ifile]]).st_size
                                if S > 0:
                                    samples[i[idx[0]]] = [i[idx[1]], os.stat(i[idx[1]]).st_size]
                                else:
                                    errors["Empty input files"].append(i[idx[ifile]])
                            except:
                                errors["Missing input files"].append(i[idx[ifile]])
                else:
                    exit("Error: Input sample files must be '.fastq' or '.fastq.gz'")
    if len(samples) == 0:
        exit("Error: No available samples identified.")
    r = 0
    for i,j in errors.iteritems():
        if len(j) > 0:
            r += 1
            print i + ":"
            for k in j:
                print "- " + k
    if r > 0:
        exit("Samples file errors detected")
    return samples


#############################################################################
## CHECKS ARGUMENTS
#############################################################################
def check_args(path_base, folder):
    if path_base == "":
        # If not provided, default path_base set to the current working directory
        path_base = os.getcwd()
    if not path_base.endswith("/"):
        path_base += "/"
    # Check if the project folder exists
    temp = os.listdir(path_base)
    folder = folder.replace("/", "")
    if not (folder in temp):
        exit("Error: Project folder '" + folder + "' not found at path_base: '" + path_base + "'")
    return [path_base, folder]


#############################################################################
# CHECKS OUTPUT FILES RELATED TO EACH ANALYSIS
#############################################################################
def check_samples(samples, path_base, folder, pg, m):
    ld = os.listdir(path_base+"/"+folder)
    if m == "update":
        if "results_"+pg in ld:
            v_samples = dict()
            ld = os.listdir(path_base+"/"+folder+"/results_"+pg)
            sok = dict()
            sko = dict()
            if os.path.exists(path_base+"/"+folder+"/results_"+pg + "/samples_ok.txt"):
                f = open(path_base+"/"+folder+"/results_"+pg + "/samples_ok.txt", 'r')
                for i in f:
                    i = i.strip("\n")
                    if len(i) > 0:
                        sok[i] = 1
                f.close()
            if os.path.exists(path_base+"/"+folder+"/results_"+pg + "/samples_ko.txt"):
                f = open(path_base+"/"+folder+"/results_"+pg + "/samples_ko.txt", 'r')
                for i in f:
                    i = i.strip("\n")
                    if len(i) > 0:
                        sko[i] = 1
                        if sok.has_key(i):
                            sok.pop(i)
                f.close()
            k = [len(samples), len(sok), len(sko)] # total, ok, ko
            for sample_name in samples.keys():
                if (not sample_name in sok.keys()) or (sample_name in sko.keys()):
                    v_samples[sample_name] = samples[sample_name]
            return v_samples, k
        else:
            return samples, [len(samples), 0, 0]
    elif m == "new":
        return samples, [len(samples), 0, 0]
    else:
        exit("Error: Mode not valid ('update' or 'new')")
