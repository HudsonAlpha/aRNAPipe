# -*- coding: utf-8 -*-
import os
import config
import vcrparser
# Conditional library load according to 'config.mode'
if config.mode == "LSF":
    import sys_LSF as manager
elif config.mode == "LOCAL":
    import sys_single as manager
else:
    import sys_OTHER as manager

sample_checker = " && (echo '#SAMPLE' >> #FOLDER/samples_ok.txt) || (echo '#SAMPLE' >> #FOLDER/samples_ko.txt)"


def submit_job_super(pname, path, wt, nproc, q, ns, bsub_suffix, nstar, timestamp):
    folder = path.split("/")[-1]
    if pname.startswith("htseq"):
        jid = "hts" + (pname.split("-")[1][0]).upper()
    elif pname.startswith("star-fusion"):
        jid = "sFusion"
    else:
        jid = pname[0:3]
    print "> Submitting " + pname + " job(s) to cluster..."
    uds = list()
    logs = list()
    for i in range(nstar):
        logname = path + "/logs/" + timestamp + "_" + pname + "_" + str(i) + ".log"
        logs.append(logname)
        if os.path.exists(logname):
            os.remove(logname)
        uds.append(str(manager.submit_job(wt, str(nproc), q, logname, folder + "_" + jid + "_" + str(i),
                              path + "/results_" + pname + "/script_" + str(i) + ".sh", 1, path, bsub_suffix)))
    uds = "|".join(uds)
    logs = "|".join(logs)
    out = open(path + "/temp/pids.txt", 'a')
    print >> out, pname + "\t" + uds + "\t" + logs
    out.close()
    return uds, logs


def secure_mkdir(path,folder):
    L = os.listdir(path)
    if not (folder in L):
        os.mkdir(path + "/" + folder)
        print "> '" + folder + "' folder created"
    else:
        print "> - '" + folder + "' folder already exists"
    return 1


def create_scripts(nchild, commands, path_base, folder, output):
    out = list()
    for i in range(nchild):
        out.append(open(path_base + folder + "/" + output + "/script_" + str(i) + ".sh",'w'))
    for k in range(len(commands)):
        print >> out[k % nchild], commands[k]
    for i in out:
        i.close()


def sortbysize(samples):
    sizes = dict()
    for sample, files in sorted(samples.iteritems()):
        if len(files) == 2:
            r = files[1]
        else:
            r = files[2] + files[3]
        if not sizes.has_key(r):
            sizes[r] = [sample]
        else:
            sizes[r].append(sample)
    keys = list()
    for i,j in sorted(sizes.iteritems()):
        for k in j:
            keys.append(k)
    return keys


def compute_mean_std(path_base, folder, samples, output, nproc, wt, q):
    ########################################################################
    ## Computes mean fragment length for kallisto single-end reads
    ########################################################################
    f  = list()
    s  = list()
    for i,j in sorted(samples.iteritems()):
        f.append("'"+j[0]+"'")
        s.append("'"+i+"'")
    s  = "samples = ["+",".join(s)+"]"
    f  = "files   = ["+",".join(f)+"]"
    o  = "output  = '"+output+"'"
    a1 = 'out = open(output, "w")\nprint >> out, "sample N M V"\nfor isamp in range(len(samples)):\n    f  = open(files[isamp],"r")\n    k  = 0'
    a2 = '    kt = 0\n    r  = list()\n    for i in f:\n        if (k % 4)==1:\n            l = float(len(i.rstrip()))\n            r.append(l)\n            kt += l\n        k += 1'
    a3 = '    f.close()\n    N = len(r)\n    M = float(kt)/N\n    t = 0\n    for i in r:\n        t += ((i-M)*(i-M))\n    S = (t/N)**0.5\n    if S < 0.001:\n        S = 0.001\n    print >> out, samples[isamp]+" "+str(int(N))+" "+str(round(M,3))+" "+str(round(S,3))'
    a4 = 'out.close()\n'
    out = open(path_base + folder + "/results_kallisto/script_stats.py",'w')
    print >> out, s
    print >> out, f
    print >> out, o
    print >> out, a1
    print >> out, a2
    print >> out, a3
    print >> out, a4
    out.close()
    print "> Submitting Kallisto stat reads job to cluster..."
    uds = manager.submit_job(wt, str(min(int(nproc),len(samples))), q, path_base+folder+"/results_kallisto/log_cluster_stats.txt", "Jmean", "python "+path_base+folder+"/results_kallisto/script_stats.py", 0, path_base+folder, "")
    print "> Kallisto stats job ID: "+uds
    return uds, path_base+folder+"/results_kallisto/log_cluster_stats.txt"


def rename_tg_output(sample, files, path):
    g = path + "/results_trimgalore/"
    cmds = list()
    if len(files) == 4:
        for i in range(2):
            output = g + files[i].split("/")[-1].replace(".fastq.gz","").replace(".fastq","") + "_val_" + str(i+1) +".fq.gz"
            cmds.append("mv " + output + " " + output.replace("_val_" + str(i+1) +".fq.gz", ".fastq.gz"))
    else:
        output = g + files[0].split("/")[-1].replace(".fastq.gz","").replace(".fastq","") + "_trimmed.fq.gz"
        cmds.append("mv " + output + " " + output.replace("_trimmed.fq.gz", ".fastq.gz"))
    return "\n".join(cmds)


def trimgalore(timestamp, path_base, folder, samples, nproc, wt, q, extra_args):
    ########################################################################
    ## FastQC analysis
    ########################################################################
    print "## Trim-galore: Quality and adapter trimming"
    print "> Quality and adapter trimming with Trim Galore..."
    output = "results_trimgalore"
    secure_mkdir(path_base + folder, "results_trimgalore")
    output_folder = path_base + folder + "/results_trimgalore"

    print "> Writing jobs for TrimGalore analysis..."
    nproc, nchild, bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    for sample in ksamp:
        files = samples[sample]
        if len(files) == 4:
            args = extra_args + " --paired"
            fnames = files[0] + " " + files[1]
        else:
            args = extra_args
            fnames = files[0]
        if (args != "") and (not args.startswith(" ")):
            args = " " + args
        call = config.path_trimgalore + args + " --gzip --path_to_cutadapt " + config.path_cutadapt + " -o " + output_folder + " " + fnames
        call = call + sample_checker.replace("#FOLDER", output_folder).replace("#SAMPLE", sample) + "\n" + rename_tg_output(sample, files, path_base + folder)
        commands.append(call)
    create_scripts(nchild, commands, path_base, folder, output)
    return submit_job_super("trimgalore", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, nchild, timestamp)


def fastqc(timestamp, path_base, folder, samples, nproc, wt, q, tg):
    ########################################################################
    ## FastQC analysis
    ########################################################################
    print "## QC: FastQC"
    print "> Quality control with fastQC..."
    output = "results_fastqc"
    secure_mkdir(path_base + folder, "results_fastqc")
    output_folder = path_base + folder + "/results_fastqc"
    print "> Writing jobs for fastqc analysis..."
    nproc, nchild, bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    for sample in ksamp:
        files = samples[sample]
        if not tg:
            if len(files) == 4:
                fnames = files[0] + " " + files[1]
            else:
                fnames = files[0]
        else:
            g = path_base + folder + "/results_trimgalore/"
            suf = ""
            if not files[0].split("/")[-1].endswith(".gz"):
                suf = ".gz"
            if len(files) == 4:
                fnames = g + files[0].split("/")[-1] + suf + " " + g + files[1].split("/")[-1] + suf
            else:
                fnames = g + files[0].split("/")[-1] + suf
        call = config.path_fastqc + " -q -o " + output_folder + " " + fnames
        commands.append(call + sample_checker.replace("#FOLDER", output_folder).replace("#SAMPLE", sample))
    create_scripts(nchild, commands, path_base, folder, output)
    return submit_job_super("fastqc", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, nchild, timestamp)


def kallisto(timestamp, path_base, folder, samples, path_index, bootstrap, nproc, wt, q, tg):
    output = "results_kallisto"
    secure_mkdir(path_base + folder, "results_kallisto")
    print "## RNAseq pseudoalignment with Kallisto"
    # Estimate counts in single-end datasss
    if len(samples[samples.keys()[0]]) == 2:
        print "> Estimating average and STD of fragment lengh required by Kalisto on single-read data..."
        outputT = path_base + folder + "/" + output + "/stats.txt"
        tid,log = compute_mean_std(path_base, folder, samples, outputT, "1", wt, q)
        vcrparser.job_wait(log, 10)
        f = open(outputT,'r')
        i = f.readline()
        stats = dict()
        for i in f:
            i = i.strip("\n").split(" ")
            stats[i[0]] = [i[2],i[3]]
        f.close()
    print "> Writing jobs for Kallisto pseudoalignment"
    nproc, nchild, bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    for sample in ksamp:
        files = samples[sample]
        if not tg:
            if len(files) == 4:
                args = ""
                fnames = files[0]+" "+files[1]
            else:
                args = " --single -l mean -s var".replace("mean", stats[sample][0]).replace("var", stats[sample][1])
                fnames = files[0]
        else:
            g = path_base + folder + "/results_trimgalore/"
            suf = ""
            if not files[0].split("/")[-1].endswith(".gz"):
                suf = ".gz"
            if len(files) == 4:
                args = ""
                fnames = g + files[0].split("/")[-1] + suf + " " + g + files[1].split("/")[-1] + suf
            else:
                args = " --single -l mean -s var".replace("mean", stats[sample][0]).replace("var", stats[sample][1])
                fnames = g + files[0].split("/")[-1] + suf
        cmd = config.path_kallisto+" quant -b " + bootstrap + " -i " + path_index + " -o " + path_base+folder + "/results_kallisto/" + sample + args + " " + fnames
        commands.append(cmd + sample_checker.replace("#FOLDER", path_base + folder + "/" + output).replace("#SAMPLE", sample))
    create_scripts(nchild, commands, path_base, folder, output)
    return submit_job_super("kallisto", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, nchild, timestamp)


def star(timestamp, path_base, folder, samples, nproc, wt, q, path_genome, star_params, tg):
    output = "results_star"
    secure_mkdir(path_base + folder, output)
    print "## RNAseq alignment with STAR..."
    print "> Writing jobs for STAR alignment..."
    nproc, nchild, bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    for sample in ksamp:
        gg = ""
        files = samples[sample]
        if not tg:
            if len(files) == 2:
                fn = files[0]
            else:
                fn = files[0] + " " + files[1]
            if files[0].endswith(".fastq.gz"):
                gg = " --readFilesCommand zcat"
        else:
            gg = " --readFilesCommand zcat"
            g = path_base + folder + "/results_trimgalore/"
            suf = ""
            if not files[0].split("/")[-1].endswith(".gz"):
                suf = ".gz"
            if len(files) == 2:
                fn = g + files[0].split("/")[-1] + suf
            else:
                fn = g + files[0].split("/")[-1] + suf + " " + g + files[1].split("/")[-1] + suf
        command = config.path_star + " --quantMode TranscriptomeSAM GeneCounts --runThreadN " + str(nproc) + " --genomeDir " + path_genome
        command = command + " --readFilesIn " + fn + " --outFileNamePrefix " + path_base + folder + "/results_star/" + sample + "_" + gg
        if len(star_params) > 0:
            command = command + star_params
        commands.append(command + sample_checker.replace("#FOLDER", path_base + folder + "/results_star").replace("#SAMPLE", sample))
    create_scripts(nchild, commands, path_base, folder, output)
    return submit_job_super("star", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, nchild, timestamp)


def starfusion(timestamp, path_base, folder, samples, nproc, wt, q, genomebuild):
    output = "results_star-fusion"
    secure_mkdir(path_base + folder, output)
    print "## Identification of gene fusions with star-fusion"
    print "> Writing jobs for Star-Fusion..."
    nproc, nchild, bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    proc_files = os.listdir(path_base + folder + "/results_star/")
    ref_file =  config.path_annotation.replace("#LABEL", genomebuild)
    for sample in ksamp:
        in_file1 = path_base + folder + "/results_star/" + sample + "_Chimeric.out.junction"
        in_file2 = path_base + folder + "/results_star/" + sample + "_Chimeric.out.sam"
        prefix = path_base + folder + "/results_star-fusion/" + sample
        if os.path.exists(in_file1) and os.path.exists(in_file2):
            call = config.path_starfusion + " -J " + in_file1 + " -S " + in_file2 + " -G " + ref_file + " --out_prefix " + prefix
            commands.append(call + sample_checker.replace("#FOLDER", path_base + folder + "/results_star-fusion").replace("#SAMPLE", sample))
        else:
            print "Warning: [Star-Fusion] STAR output file not found -> " + in_file1
    create_scripts(nchild, commands, path_base, folder, output)
    return  submit_job_super("star-fusion", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, nchild, timestamp)


def picardqc(timestamp, path_base, folder, samples, nproc, wt, q, annots, strand):
    nstrand = {" --stranded=no":"NONE", " --stranded=yes":"FIRST_READ_TRANSCRIPTION_STRAND", " --stranded=no":"SECOND_READ_TRANSCRIPTION_STRAND"}
    output = "results_picard"
    secure_mkdir(path_base + folder, output)
    print "## Alignment QC Picard-CollectRnaSeqMetrics"
    print "> Writing jobs for Picard QC..."
    nproc, nchild, bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    proc_files = os.listdir(path_base + folder + "/results_star/")
    for sample in ksamp:
        in_file = path_base + folder + "/results_star/" + sample + "_Aligned.out.sam"
        if sample + "_Aligned.out.sam" in proc_files:
            for i in range(len(config.nannots)):
                annot    = annots[i]
                out_file = in_file.replace(".sam", "." + config.nannots[i] + ".qc").replace("results_star/", "results_picard/").replace("_Aligned.out", "")
                call = "java -jar " + config.path_picard + "/CollectRnaSeqMetrics.jar REF_FLAT=" + annot + " STRAND_SPECIFICITY=" + nstrand[strand] + " INPUT=" + in_file + " OUTPUT=" + out_file
                if i == (len(config.nannots)-1):
                    commands.append(call + sample_checker.replace("#FOLDER", path_base + folder + "/results_picard").replace("#SAMPLE", sample))
                else:
                    commands.append(call)
        else:
            print "Warning: [Picard] STAR output file not found -> " + in_file
    create_scripts(nchild, commands, path_base, folder, output)
    return  submit_job_super("picard", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, nchild, timestamp)


def htseq(timestamp, path_base, folder, samples, path_annotation, nproc, wt, q, mode, strand, countmode):
    output = "results_htseq-" + mode
    secure_mkdir(path_base + folder, output)
    print "## HTseq-count"
    print "> Writing jobs for HTseq-count " + mode + " analysis..."
    nproc, nchild, bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    proc_files = os.listdir(path_base + folder + "/results_star/")
    for sample in ksamp:
        in_file = path_base + folder + "/results_star/" + sample + "_Aligned.out.sam"
        if sample + "_Aligned.out.sam" in proc_files:
            outputf= path_base + folder + "/results_htseq-" + mode + "/" + sample + ".tab"
            if mode == "gene":
                ld1 = config.path_htseq + strand + " -m " + countmode  + " -q " + in_file + " " + path_annotation
            else:
                ld1 = config.path_htseq + strand + " -m " + countmode  + " -i exon_id -q " + in_file + " " + path_annotation
            call = ld1 + " > " + outputf
            commands.append(call  + sample_checker.replace("#FOLDER", path_base + folder + "/" + output).replace("#SAMPLE", sample))
        else:
            print "Warning: [HTseq-" + mode + "] STAR output file not found -> " + in_file
    create_scripts(nchild, commands, path_base, folder, output)
    return  submit_job_super("htseq-" + mode, path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, nchild, timestamp)


def sam2sortbam(timestamp, path_base, folder, samples, nproc, wt, q):
    output = "results_sam2sortbam"
    secure_mkdir(path_base + folder, output)
    print "## SAM2SORTEDBAM"
    print "> Writing jobs for SAM2SORTEDBAM..."
    nproc, nchild, bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    proc_files = os.listdir(path_base + folder + "/results_star/")
    for sample in ksamp:
        in_file = path_base + folder + "/results_star/" + sample + "_Aligned.out.sam"
        if sample + "_Aligned.out.sam" in proc_files:
            out_file = path_base + folder + "/results_sam2sortbam/" + sample + ".sorted.bam"
            com = "java -jar " + config.path_picard + "/AddOrReplaceReadGroups.jar I=" + in_file + " O=" + out_file +" SO=coordinate RGID=id RGLB=library RGPL=ILLUMINA RGPU=machine RGSM=sample 2> " + out_file + ".log"
            commands.append(com + sample_checker.replace("#FOLDER", path_base + folder + "/results_sam2sortbam").replace("#SAMPLE", sample))
        else:
            print "Warning: [SAM2SORTEDBAM] STAR output file not found -> " + in_file
    create_scripts(nchild, commands, path_base, folder, output)
    return  submit_job_super("sam2sortbam", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, nchild, timestamp)


def jsplice(timestamp, path_base, folder, samples, nproc, wt, q, genomebuild, pheno):
    output_dir = path_base + folder + '/results_jsplice'
    secure_mkdir(path_base + folder, output_dir)
    print "## jSPLICE"
    print "> Writing jobs for jSPLICE..."
    nproc, nchild, bsub_suffix = manager.get_bsub_arg('1/NA/NA', len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    proc_files_1 = os.listdir(path_base + folder + "/results_star/")
    proc_files_2 = os.listdir(path_base + folder + "/results_sam2sortbam/")
    out = open(output_dir + '/expdesign.txt', 'w')
    print >> out, '#exp\tcond\tjxnFile\tbamFile'
    for sample in ksamp:
        sj_file = path_base + folder + '/results_star/' + sample + '_SJ.bed' # Junction file created by STAR
        sj_out_file = output_dir + '/' + sample + '.SJ.bed'
        bam_file = path_base + folder + '/results_sam2sortbam/' + sample + '.sorted.bam' # BAM file created by STAR/Picard(AddOrReplaceReadGroups)
        if (sj_file in proc_files_1) and (bam_file in proc_files_2):
            command = 'python ' + config.path_jsplice + '/starJxn2bed.py -f ' + sj_file + ' -o '+ sj_out_file
            commands.append(command + sample_checker.replace("#FOLDER", output_dir).replace("#SAMPLE", sample))
            print >> out, '\t'.join([pheno[sample].split(':')[0], pheno[sample].split(':')[1], sj_out_file, bam_file])
        else:
            print "Warning: [JSPLICE] STAR output files not found -> " + sample
    out.close()
    commands.append('python ' + config.path_jsplice + '/jSplice.py  –d ' + output_dir + '/expdesign.txt –o ' + output_dir + ' -a '+ config.path_annotation.replace("#LABEL", genomebuild) + ' -c 10')
    create_scripts(nchild, commands, path_base, folder, 'results_jsplice')
    return  submit_job_super("sam2sortbam", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, nchild, timestamp)


def picard_IS(timestamp, path_base, folder, samples, nproc, wt, q):
    output = "results_picard_IS"
    secure_mkdir(path_base + folder, output)
    print "## Picard-InsertSize"
    print "> Writing jobs for Picard InsertSize..."
    nproc, nchild, bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    proc_files = os.listdir(path_base + folder + "/results_sam2sortbam/")
    for sample in ksamp:
        in_file = path_base + folder + "/results_sam2sortbam/" + sample + ".sorted.bam"
        if sample + ".sorted.bam" in proc_files:
            for i in range(len(config.nannots)):
                out_file = in_file.replace("results_sam2sortbam/", "results_picard_IS/").replace(".sorted.bam", "")
                call = "java -jar " + config.path_picard + "/CollectInsertSizeMetrics.jar I="+in_file+" O="+out_file+".txt H="+out_file+".pdf"
                commands.append(call + sample_checker.replace("#FOLDER", path_base + folder + "/results_picard_IS").replace("#SAMPLE", sample))
        else:
            print "Warning: [Picard] Sorted BAM file not found -> " + in_file
    create_scripts(nchild, commands, path_base, folder, output)
    return  submit_job_super("picard_IS", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, nchild, timestamp)


def varscan(timestamp, path_base, folder, samples, nproc, wt, q, genome_build, args):
    ref = config.path_fasta.replace("#LABEL",genome_build)
    output = "results_varscan"
    secure_mkdir(path_base + folder, output)
    print "## Variang calling with VARSCAN"
    print "> Writing jobs for VARSCAN..."
    nproc, nchild, bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    proc_files = os.listdir(path_base + folder + "/results_sam2sortbam/")
    for sample in ksamp:
        in_file = path_base + folder + "/results_sam2sortbam/" + sample + ".sorted.bam"
        if sample + ".sorted.bam" in proc_files:
            out_file = path_base + folder + "/results_varscan/" + sample + ".vcf"
            com = config.path_samtools + " mpileup -B -f " + ref + " " + in_file + " | java -jar " + config.path_varscan + " mpileup2cns " + args + " > " + out_file
            commands.append(com + sample_checker.replace("#FOLDER", path_base + folder + "/results_varscan").replace("#SAMPLE", sample))
        else:
            print "Warning: [VARSCAN] SORTED BAM output file not found -> " + in_file
    create_scripts(nchild, commands, path_base, folder, output)
    return  submit_job_super("varscan", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, nchild, timestamp)


def gatk(timestamp, path_base, folder, samples, nproc, wt, q, genome_build, args):
    args = args.split("|")
    multithread = False
    filt = "30"
    if len(args) == 2:
        if args[0] == "yes":
            multithread = True
        filt = args[1]
    output = "results_gatk"
    secure_mkdir(path_base + folder, output)
    print "## Variang calling with GATK"
    print "> Writing jobs for GATK..."
    nproc, nchild, bsub_suffix = manager.get_bsub_arg(nproc, len(samples))
    commands = list()
    ksamp = sortbysize(samples)
    proc_files = os.listdir(path_base + folder + "/results_sam2sortbam/")
    for sample in ksamp:
        in_file = path_base + folder + "/results_sam2sortbam/" + sample + ".sorted.bam"
        if sample + ".sorted.bam" in proc_files:
            C = gatk_commands(path_base + folder, sample, genome_build, multithread, filt)
            commands.append("\n".join(C))
        else:
            print "Warning: [GATK] SORTED BAM output file not found -> " + in_file
    create_scripts(nchild, commands, path_base, folder, output)
    return  submit_job_super("gatk", path_base + folder, wt, str(nproc), q, len(samples), bsub_suffix, nchild, timestamp)


def gatk_commands(path, sample, ref, multithread, filter):
    try:
        parameters = { "SAMPLE": sample,
                       "INPUT": path + "/results_sam2sortbam/" + sample + ".sorted.bam",
                       "OUTDIR": path + "/results_gatk/",
                       "PICARD": config.path_picard,
                       "GATK": config.path_gatk,
                       "MT_RTC": "",
                       "MT_BR": "",
                       "MT_PR": "",
                       "FILTERVAL": str(filter),
                       "REF": config.path_fasta.replace("#LABEL", ref),
                       "DBSNP": config.path_db + "/genomes_processed/" + ref + "/gatk/" + config.annots_gatk[ref][0]}
        temp = ["", ""]
        for i in config.annots_gatk[ref][1]:
            temp[0] = temp[0] + " -known " + config.path_db + "/genomes_processed/" + ref + "/gatk/" + i
            temp[1] = temp[1] + " -knownSites " + config.path_db + "/genomes_processed/" + ref + "/gatk/" + i
        parameters["INDEL"] = temp[0]
        parameters["IND_SITES"] = temp[1]
        if multithread:
            parameters["MT_RTC"] = " -nt " + str(config.gatk_multithread["RTC"])
            parameters["MT_BR"] = " -nct " + str(config.gatk_multithread["BR"])
            parameters["MT_PR"] = " -nct " + str(config.gatk_multithread["PR"])
    except:
        return []
    C = list()
    # MARK DUPLICATES AND CREATE INDEX
    C.append("java -jar #PICARD/MarkDuplicates.jar  I=#INPUT O=#OUTDIR/#SAMPLE.dedup.bam CREATE_INDEX=true VALIDATION_STRINGENCY=SILENT M=#OUTDIR/#SAMPLE.sorted.dedup.metrics 2> #OUTDIR/#SAMPLE.1.dup.log")
    # SPLIT N CIGARS
    C.append("java -jar #GATK -T SplitNCigarReads -R #REF -I #OUTDIR/#SAMPLE.dedup.bam -o #OUTDIR/#SAMPLE.split.bam -rf ReassignOneMappingQuality -RMQF 255 -RMQT 60 -U ALLOW_N_CIGAR_READS 2> #OUTDIR/#SAMPLE.2.split.log")
    # INDEL REALIGNMENT
    C.append("java -jar #GATK -T RealignerTargetCreator #MT_RTC -R #REF -I #OUTDIR/#SAMPLE.split.bam -o #OUTDIR/#SAMPLE.intervals #INDEL 2> #OUTDIR/#SAMPLE.3.realign1.log")
    C.append("java -jar #GATK -T IndelRealigner -R #REF -I #OUTDIR/#SAMPLE.split.bam -targetIntervals #OUTDIR/#SAMPLE.intervals #INDEL -o #OUTDIR/#SAMPLE.processed.bam 2> #OUTDIR/#SAMPLE.4.realign2.log")
    # BASE RECALIBRATION
    C.append("java -jar #GATK -T BaseRecalibrator #MT_BR -I #OUTDIR/#SAMPLE.processed.bam -R #REF #IND_SITES -knownSites #DBSNP -o #OUTDIR/#SAMPLE.table 2> #OUTDIR/#SAMPLE.5.BQSR.log")
    C.append("java -jar #GATK -T PrintReads #MT_PR -R #REF -I #OUTDIR/#SAMPLE.processed.bam -BQSR #OUTDIR/#SAMPLE.table -o #OUTDIR/#SAMPLE.rec.bam 2> #OUTDIR/#SAMPLE.6.BQSR2.log")
    # HAPLOTYPE CALLER
    C.append("java -jar #GATK -T HaplotypeCaller -R #REF -I #OUTDIR/#SAMPLE.rec.bam -dontUseSoftClippedBases -stand_call_conf 20.0 -stand_emit_conf 20.0 -o #OUTDIR/#SAMPLE.vcf 2> #OUTDIR/#SAMPLE.7.HaploCall.log")
    # FILTER
    C.append("java -jar #GATK -T VariantFiltration -R #REF -V #OUTDIR/#SAMPLE.vcf -window 35 -cluster 3 -filterName FS -filter 'FS > #FILTERVAL' -filterName QD -filter 'QD < 2.0' -o #OUTDIR/#SAMPLE.filt.vcf 2>  #OUTDIR/#SAMPLE.8.filt.log")
    # SUBSTITUTION
    for label, value in parameters.iteritems():
        for i in range(len(C)):
            C[i] = C[i].replace("#" + label, value)
    C[len(C)-1] = C[len(C)-1] + sample_checker.replace("#FOLDER", path + "/results_gatk").replace("#SAMPLE", sample)
    return C

