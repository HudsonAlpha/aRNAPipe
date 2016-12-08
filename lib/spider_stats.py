import os
import config

head_star_log = ["Sample", "Links","Started job","Started mapping","Finished","Mapping speed [Mr/h]","Input reads [n]","Input read length (mean)","Uniquely mapped [n]","Uniquely mapped [%]","Mapped length","Splices [n]","Splices annotated [n]","Splices GT/AG [n]","Splices: GC/AG [n]","Splices: AT/AC [n]","Splices: Non-canonical [n]","Mismatch rate per base [%]","Deletion rate per base [%]","Deletion average length","Insertion rate per base [%]","Insertion average length","Multimapping reads [n]","Multimapping reads [%]","Multimapping reads (+) [n]","Multimapping reads (+) [%]","Unmapped reads: too many mismatches [%]","Unmapped reads: too short [%]","Unmapped reads: other [%]"]


def stats_varcall(path, mod):
    f = open(path + "/samples.list", 'r')
    h = f.readline()
    samples = dict()
    for i in f:
        i = i.strip("\n").split("\t")
        samples[i[0]] = path + "/results_" + mod + "/" + i[0] + ".vcf"
    f.close()
    out = open(path + "/outputs/stats_" + mod + ".txt", 'w')
    print >> out, "sample_id\t0/1\t1/1"
    for sample, filename in samples.iteritems():
        if os.path.exists(filename):
            f = open(filename, 'r')
            i = f.readline()
            n = [0,0,0]
            for i in f:
                if not i.startswith("#"):
                    i = i.strip("\n").split("\t")
                    if len(i) > 1:
                        gt = i[9].split(":")[0]
                        if gt == "1/1":
                            n[1] += 1
                        elif gt == "0/1":
                            n[0] += 1
                        elif gt == "1/2":
                            n[2] += 1
                        else:
                            print gt
            f.close()
            print >> out, sample + "\t" + str(n[0]) + "\t" + str(n[1]) + "\t" + str(n[2])
        else:
            print >> out, sample + "\t-1\t-1\t-1"
    out.close()

def stats_fusion(path):
    f = open(path + "/samples.list", 'r')
    h = f.readline()
    samples = dict()
    for i in f:
        i = i.strip("\n").split("\t")
        samples[i[0]] = path + "/results_star-fusion/" + i[0] + "/star-fusion.fusion_candidates.final.abridged"
    f.close()
    out = open(path + "/outputs/starfusion_aggregate.txt", 'w')
    print >> out, "sample_id\tFusionName\tJunctionReadCount\tSpanningFragCount\tSpliceType\tLeftGene\tLeftBreakpoint\tRightGene\tRightBreakpoint\tLargeAnchorSupport\tLeftBreakDinuc\tLeftBreakEntropy\tRightBreakDinucRightBreakEntropy"
    for sample, filename in samples.iteritems():
        if os.path.exists(filename):
            f = open(filename, 'r')
            i = f.readline()
            for i in f:
                i = i.strip("\n").split("\t")
                if len(i) > 1:
                    print >> out, sample + "\t" + "\t".join(i)
    out.close()


def stats_trimgalore(path):
    fields2 = ["Total reads processed:","Reads with adapters:","Reads written (passing filters):","Total basepairs processed:","Quality-trimmed:","Total written (filtered):"]
    fnames2 = ['Processed reads', 'Reads with adapters', 'Reads passing filters', 'Processed basepairs', 'Quality-trimmed basepairs', 'Basepairs passing filters']
    f = open(path + "/samples.list", 'r')
    h = f.readline().strip('\n').split('\t')
    idx = []
    for ix, i in enumerate(h):
        if i.startswith('FASTQ'):
            idx.append(ix)
    samples = dict()
    for i in f:
        i = i.strip("\n").split("\t")
        if len(i) > 1:
            samples[i[0]] = [j.split("/")[-1] for j in idx]
    f.close()
    out = open(path + "/outputs/stats_trim.txt", 'w')
    print >> out, "sample_id\t" + "\t".join(fnames2)
    out2 = open(path + "/outputs/stats_trim_plot.txt", 'w')
    print >> out2, "sample_id\tReads with adapters\tReads passing filters\tQuality-trimmed basepairs\tBasepairs passing filters"
    for sample, files in samples.iteritems():
        k = 0
        data = [{},{}]
        for filname in files:
            fpath = path + "/results_trimgalore/" + filname + "_trimming_report.txt"
            if os.path.exists(fpath):
                print 'open'
                f = open(fpath, 'r')
                for i in f:
                    i = i.strip("\n").replace(' bp', '').replace(',', '')
                    for fi in fields2:
                        j = i.split(fi)
                        if len(j) > 1:
                            data[k][fi] = ' '.join(j[1:])
                f.close()
            k += 1
        g1 = ""
        g2 = ""
        plt1 = ''
        plt2 = ''
        for i in fields2:
            r1 = "NA"
            r2 = "NA"
            if data[0].has_key(i):
                r1 = data[0][i]
            if len(files) == 2:
                if data[1].has_key(i):
                    r2 = data[1][i]
                g1 = g1 + "\t" + r1
                g2 = g2 + "\t" + r2
                if '%' in r2:
                    plt1 = plt1 + "\t" + r1.split('(')[1].replace('%)', '')
                    plt2 = plt2 + "\t" + r2.split('(')[1].replace('%)', '')
            else:
                if '%' in r1:
                    plt1 = plt1 + "\t" + r1.split('(')[1].replace('%)', '')
                g1 = g1 + "\t" + r1
        if g2 != '':
            print >> out, sample + ' (1)' + g1
            print >> out, sample + ' (2)' + g2
            print >> out2, sample + ' (1)' + plt1
            print >> out2, sample + ' (2)' + plt2
        else:
            print >> out, sample + g1
            print >> out2, sample + plt1
    out.close()


def read_star_log(datafile, sample):
    H = []
    if os.path.exists(datafile):
        link1 = "../results_star/" + sample + "_Log.final.out"
        link1 = '<a href="LINK" target="_blank">FINAL</a>'.replace("LINK",link1)
        link2 = "../results_star/"+sample+"_Log.progress.out"
        link2 = '<a href="LINK" target="_blank">PROGRESS</a>'.replace("LINK",link2)
        link3 = "../results_star/"+sample+"_ReadsPerGene.out.tab"
        link3 = '<a href="LINK" target="_blank">COUNTS</a>'.replace("LINK",link2)
        links = link1 + " / " + link2 + " / " + link3
        ff = open(datafile, 'r')
        n = list()
        for i in ff:
            if "|" in i:
                i = i.strip("\n").split()
                if "%" in i[len(i)-1]:
                    i[len(i)-1] = i[len(i)-1].replace("%", "")
                n.append(i[len(i)-1])
        ff.close()
        return [sample, links] + n
    else:
        n = [-1 for i in range(len(head_star_log))]
        return [sample, "-1"] + n


# def get_star_log(path, samples):
#     print "> Recovering stats from STAR logs: " + path
#     g = os.listdir(path)
#     out = open(path.replace("/results_star/", "/outputs/star_stats_log.txt") , "w")
#     print >> out, "\t".join(head_star_log)
#     for sample in samples:
#         if os.path.exists(path + "/" + sample + "_Log.final.out"):
#             n = read_star_log(path + "/" + sample + "_Log.final.out", sample)
#             print >> out, "\t".join(n)
#     out.close()



def stats_fastq(path,samples,config):
    if not os.path.exists(path):
        return 1
    n = os.listdir(path)
    if config.has_key("fastqc") and ("results_fastqc" in n):
        files = os.listdir(path+"/results_fastqc")
        h = ["Sample",
             "Link",
             '<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/1%20Basic%20Statistics.html" target="_blank">Basic statistics</a>',
             '<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/2%20Per%20Base%20Sequence%20Quality.html" target="_blank">Per base sequence quality</a>',
             '<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/3%20Per%20Sequence%20Quality%20Scores.html" target="_blank">Per sequence quality scores</a>',
             '<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/4%20Per%20Base%20Sequence%20Content.html" target="_blank">Per base sequence content</a>',
             '<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html" target="_blank">Per base GC content</a>',
             '<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/5%20Per%20Sequence%20GC%20Content.html" target="_blank">Per sequence GC content</a>',
             '<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/6%20Per%20Base%20N%20Content.html" target="_blank">Per base N content</a>',
             '<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/7%20Sequence%20Length%20Distribution.html" target="_blank">Sequence Length Distribution</a>',
             '<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/8%20Duplicate%20Sequences.html" target="_blank">Sequence Duplication Levels</a>',
             '<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/9%20Overrepresented%20Sequences.html" target="_blank">Overrepresented sequences</a>',
             '<a href="http://www.bioinformatics.babraham.ac.uk/projects/fastqc/Help/3%20Analysis%20Modules/11%20Kmer%20Content.html" target="_blank">Kmer Content</a>',
             "Total sequences",
             "Sequence length",
             "%GC"]
        names = ["","","per_base_quality.png","per_sequence_quality.png","per_base_sequence_content.png","per_base_gc_content.png","per_sequence_gc_content.png",
                 "per_base_n_content.png", "sequence_length_distribution.png", "duplication_levels.png", "", "kmer_profiles.png", "", "", ""]
        table = list()
        n = ""
        for i in h:
            n = n+"<th bgcolor='#A8A8A8'>"+i+"</th>"
        n = "<tr>"+n+"</tr>"
        table.append(n)
        iin = 0
        for sample,data in sorted(samples.iteritems()):
            filess = data[0:(len(data)/2)]
            m = list()
            for f in filess:
                n = list()
                f = f.split("/")
                f = f[len(f)-1]
                if f.replace(".fastq","").replace(".gz","")+"_fastqc" in files:
                    link = "../results_fastqc/"+f.replace(".fastq","").replace(".gz","")+"_fastqc/fastqc_report.html"
                    link = '<a href="LINK" target="_blank">Results</a>'.replace("LINK",link)
                    n.append(link)
                    ff = open(path+"/results_fastqc/"+f.replace(".fastq","").replace(".gz","")+"_fastqc/summary.txt",'r')
                    for i in ff:
                        i = i.strip("\n").split("\t")
                        if len(names[len(n)]) > 0:
                            G = '<a href="../results_fastqc/'+f.replace(".fastq","").replace(".gz","")+"_fastqc"+'/Images/'+names[len(n)]+'" class="lytebox lytetip" data-tip="" data-lyte-options="showPrint:true tipStyle:classic" data-lightbox="image-'+str(iin)+'" data-title="#TIT">#VAL</a>'
                            iin += 1
                        else:
                            G = "#VAL"
                        n.append(G.replace("#VAL",i[0]).replace("#TIT", sample + "_" + str(len(m)+1)))
                    ff.close()
                    ff = open(path+"/results_fastqc/"+f.replace(".fastq","").replace(".gz","")+"_fastqc/fastqc_data.txt",'r')
                    k=0
                    for i in ff:
                        if k==6:
                            n.append(i.strip("\n").split("\t")[1])
                        if k==8:
                            n.append(i.strip("\n").split("\t")[1])
                        if k==9:
                            n.append(i.strip("\n").split("\t")[1])
                        k+=1
                        if k>10:
                            break
                    ff.close()
                else:
                    n = ["NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA","NA"]
                m.append(n)
            s = ["<td bgcolor='#A8A8A8'>"+sample+"</td>"]
            for i in range(len(m[0])):
                if len(m)>1:
                    g = m[0][i]+" / "+m[1][i]
                else:
                    g = m[0][i]
                if ("NA" in g) or ("FAIL" in g):
                    cl = "#CC3300"
                elif "WARN" in g:
                    cl = "#FFCC00"
                else:
                    cl = "#00CC66"
                s.append("<td bgcolor='"+cl+"'>"+g+"</td>")
            s = "<tr>"+"".join(s)+"</tr>"
            table.append(s)
        return "<table>" + "\n".join(table) + "</table>"
    else:
        return ""


def get_annotation(feature, path):
    aRNApipe_path = config.path_db + "/genomes_processed/"
    genome = ""
    f = open(path, 'r')
    for i in f:
        i = i.strip("\n").split("\t")
        if len(i) > 1:
            if i[0] == "genome_build":
                genome = i[1]
                break
    f.close()
    if genome == "":
        exit("Error recovering genome build")
    path = aRNApipe_path + "/" + genome + "/genesets.feature.txt".replace("feature", feature)
    if not os.path.exists(path):
        exit("Annotation file not found")
    f = open(path, 'r')
    h = f.readline().strip("\n").split("\t")
    k = 0
    for i in range(len(h)):
        if h[i] == feature:
            k = i
            break
    d = dict()
    L = dict()
    for i in f:
        i = i.strip("\n").split("\t")
        d[i[k]] = "\t".join(i)
        L[i[k]] = i[6]
    f.close()
    return d, h, L


def parselog(path):
    # Given a LSF log file returns a dictionary with the relevant parameters
    M = {"CPU time :":"cpu_time","Max Memory :":"max_memo",
             "Average Memory :":"ave_memo","Max Swap :":"max_swap",}
    f = open(path,'r')
    results = dict()
    for i in f:
        if i.startswith("Started at"):
            results["time_start"] = i.rstrip().replace("Started at ","")
        elif i.startswith("Results reported on"):
            results["time_results"] = i.rstrip().replace("Results reported on ","")
        elif i.startswith("Successfully completed."):
            results["success"] = "Yes"
        else:
            for j in M.keys():
                if j in i:
                    results[M[j]] = " ".join(i.rstrip().split()[3:5])
    if not results.has_key("success"):
        results["success"] = "No"
    path = path.split("/")
    results["link"] = '<a href="LINK" target="_blank">+</a>'.replace("LINK","../logs/" + path[-1])
    return results


def stats_log(path):
    # Given a path to a 'logs' folder parses all the LSF log files and writes the relevant parameters
    # in a table format on an ouptut file
    order = ["success", "time_start", "time_results", "cpu_time", "max_memo", "ave_memo", "max_swap", "link"]
    Lorder= ["Success", "Started", "Ended", "CPU time", "Memory (max) [MB]", "Memory (mean) [MB]", "Swap (max) [MB]", "Link"]
    progs = ["trimgalore", "fastqc", "kallisto", "star", "star-fusion", "picard",
             "htseq-gene", "htseq-exon", "sam2sortbam", "picard_IS", "varscan", "gatk", "jsplice"]
    print "> Recovering LSF stats from: " + path
    L = os.listdir(path)
    logs = dict()
    for i in L:
        if i.endswith(".log") and (not i.endswith("_M.log")):
            j = i.replace(".log","").replace('picard_IS', 'picard-IS').split("_")
            if (len(j) == 4):
                date, time, prog, proc = j
                prog = 'picard_IS' if prog == 'picard-IS' else prog
                if not logs.has_key(date + "_" + time):
                    logs[date + "_" + time] = dict()
                if not logs[date + "_" + time].has_key(prog):
                    logs[date + "_" + time][prog] = dict()
                logs[date + "_" + time][prog][int(proc)] = [path + "/" + i, parselog(path + "/" + i)]
    out = open(path.replace("/logs/", "/outputs/") + "log_stats.txt", 'w')
    print >> out, "Run number\tTimestamp\tModule\tJob\t"+"\t".join(Lorder)
    if len(logs) > 0:
        ki = 0
        k = 0
        for timestamp in sorted(logs.keys()):
            ki += 1
            for prog in progs:
                if logs[timestamp].has_key(prog):
                    for proc, data in (logs[timestamp][prog].iteritems()):
                        n = [str(ki), timestamp, prog, str(proc)]
                        for ordi in order:
                            if data[1].has_key(ordi):
                                n.append(data[1][ordi])
                            else:
                                n.append("NA")
                        print >> out, ("\t".join(n)).replace(" sec.", "").replace(" MB", "")
                        k += 1
        print "  Recovered stats from " + str(k) + " log files."
    out.close()
    return 1


def stats_kallisto(path, samples):
    print "> Recovering stats from Kallisto and building count matrices: " + path
    M = ["eff_length", "est_counts", "tpm"]
    if not os.path.exists(path):
        return 1
    g = os.listdir(path)
    k = 0
    N = [{},{},{}]
    ng = list()
    for sample in samples:
        if sample in g:
            if os.path.exists(path + "/" + sample + "/abundance.tsv"):
                for i in range(3):
                    N[i][sample] = []
                f = open(path + "/" + sample + "/abundance.tsv", 'r')
                i = f.readline()
                nt= list()
                for i in f:
                    i = i.strip("\n").split("\t")
                    nt.append(i[0] + "\t" + i[1])
                    for k in range(3):
                        N[k][sample].append(i[k+2])
                f.close()
                if len(nt) > len(ng):
                    ng = nt
    for i in range(len(M)):
        MT = dict()
        out = open(path.replace("/results_kallisto/","/outputs/") + "/kallisto_" + M[i] + ".txt", 'w')
        print >> out, "target_id\tlength\t" + "\t".join(samples)
        for j in range(len(ng)):
            r = list()
            for sample in samples:
                if not MT.has_key(sample):
                    MT[sample] = 0
                if N[i].has_key(sample):
                    if j < len(N[i][sample]):
                        r.append(N[i][sample][j])
                        MT[sample] += float(N[i][sample][j])
                    else:
                        r.append("NA")
                else:
                    r.append("NA")
            print >> out, ng[j] + "\t" + "\t".join(r)
        out.close()
        out = open(path.replace("/results_kallisto/","/outputs/") + "/kallisto_stats_" + M[i] + ".txt", 'w')
        print >> out, "sample_id\tTotalCounts\tLink"
        for sample in samples:
            if N[i].has_key(sample):
                print >> out, sample + "\t" + str(MT[sample]) + "\t" + '<a href="../results_kallisto/'+sample +'/abundance.tsv" target="_blank">+</a>'
            else:
                print >> out, sample + "\t-1\t-1"
        out.close()

    print "  Data found for " + str(len(N[0])) + " of " + str(len(samples)) + " samples"
    print "  Creating annotation file..."
    D, H, L = get_annotation("transcript", path.replace("/results_kallisto/", "/config.txt"))
    out = open(path.replace("/results_kallisto/","/outputs/") + "/kallisto_annotation.txt", 'w')
    print >> out, "\t".join(H)
    k = 0
    for transcript in ng:
        if D.has_key(transcript.split("\t")[0].split(".")[0]):
            print >> out, D[transcript.split("\t")[0].split(".")[0]]
            k += 1
        else:
            print >> out, "\t".join(["NA" for i in range(len(H))])
    out.close()
    print "  Recovered annotation for " + str(k) + " of " + str(len(ng)) + " transcripts"
    return 1


def stats_star(path, samples):
    print "> Recovering stats from STAR logs: " + path
    if not os.path.exists(path):
        return 1
    g = os.listdir(path)
    out = open(path.replace("/results_star/", "/outputs/star_stats_log.txt") , "w")
    print >> out, "\t".join(head_star_log)
    for sample in samples:
        if os.path.exists(path + "/" + sample + "_Log.final.out"):
            n = read_star_log(path + "/" + sample + "_Log.final.out", sample)
            print >> out, "\t".join(n)
    out.close()
    print "> Recovering stats from STAR counts and building count matrices: " + path
    M = ["unstranded", "stranded", "stranded-reverse"]
    g = os.listdir(path)
    k = 0
    N = [{},{},{}]
    S = [{},{},{}]
    ng = list()
    na = list()
    for sample in samples:
        if sample + "_ReadsPerGene.out.tab" in g:
            for i in range(3):
                N[i][sample] = []
                S[i][sample] = []
            f = open(path + "/" + sample  + "_ReadsPerGene.out.tab", 'r')
            # header stats
            naa = list()
            for i in range(4):
                i = f.readline().strip("\n").split("\t")
                naa.append(i[0])
                for k in range(3):
                    S[k][sample].append(i[k+1])
            if len(naa) > len(na):
                na = naa
            # counts
            nt = list()
            for i in f:
                i = i.strip("\n").split("\t")
                nt.append(i[0])
                for k in range(3):
                    N[k][sample].append(i[k+1])
            f.close()
            if len(nt) > len(ng):
                ng = nt
            if len(naa) > len(na):
                na = naa
    for i in range(len(M)):
        MT = dict()
        out = open(path.replace("/results_star/","/outputs/") + "/star_" + M[i] + "_counts.txt", 'w')
        print >> out, "gene_id\t" + "\t".join(samples)
        for j in range(len(ng)):
            r = list()
            for sample in samples:
                if not MT.has_key(sample):
                    MT[sample] = 0
                if N[i].has_key(sample):
                    if j < len(N[i][sample]):
                        r.append(N[i][sample][j])
                        MT[sample] += int(N[i][sample][j])
                    else:
                        r.append("NA")
                else:
                    r.append("NA")
            print >> out, ng[j] + "\t" + "\t".join(r)
        out.close()
        out = open(path.replace("/results_star/","/outputs/") + "/star_" + M[i] + "_stats.txt", 'w')
        print >> out, "sample_id\t" + "\t".join(na) + "\tTotalCounts\tLink"
        for sample in samples:
            if S[i].has_key(sample):
                if len(S[i][sample]) > 0:
                    print >> out, sample + "\t" + "\t".join(S[i][sample]) + "\t" + str(MT[sample]) + "\t" +  '<a href="../results_star/'+sample+'_ReadsPerGene.out.tab" target="_blank">+</a>'
                else:
                    print >> out, sample + "\t" + "\t".join(["0" for k in range(len(na))]) + "\t0\t"
            else:
                print >> out, sample + "\t" + "\t".join(["0" for k in range(len(na))]) + "\t0\t"
        out.close()
    print "  Data found for " + str(len(N[0])) + " of " + str(len(samples)) + " samples"

    print "  Creating annotation file..."
    D, H, L = get_annotation("gene", path.replace("/results_star/", "/config.txt"))
    out = open(path.replace("/results_star/","/outputs/") + "/star_annotation.txt", 'w')
    print >> out, "\t".join(H)
    k = 0
    for gene in ng:
        if D.has_key(gene):
            print >> out, D[gene]
            k += 1
        else:
            print >> out, "\t".join(["NA" for i in range(len(H))])
    out.close()
    print "  Recovered annotation for " + str(k) + " of " + str(len(ng)) + " genes"
    print "  Computing RPKMs..."
    for i in range(len(M)):
        out = open(path.replace("/results_star/","/outputs/") + "/star_" + M[i] + "_rpkms.txt", 'w')
        print >> out, "gene_id\t" + "\t".join(samples)
        for j in range(len(ng)):
            r = list()
            for sample in samples:
                if N[i].has_key(sample):
                    if (len(N[i][sample]) > j) and (MT[sample] > 0):
                        s = str(round(float(N[i][sample][j]) * 1000000000 / (float(L[ng[j]]) * MT[sample]),4))
                    else:
                        s = "NA"
                    r.append(s)
                else:
                    r.append("NA")
            print >> out, ng[j] + "\t" + "\t".join(r)
        out.close()
    return 1


def stats_htseq(path, samples, mode):
    print "> Recovering stats from HTseq " + mode + " counts and building count matrices: " + path
    special = "no_feature","ambiguous","too_low_aQual","not_aligned","alignment_not_unique"
    if not os.path.exists(path):
        return 0,0
    g = os.listdir(path)
    k = 0
    N = [{}]
    S = [{}]
    na = list()
    ng = list()
    for sample in samples:
        if sample + ".tab" in g:
            for i in range(1):
                N[i][sample] = []
                S[i][sample] = []
            f = open(path + "/" + sample  + ".tab", 'r')
            naa = list()
            ngg = list()
            for i in f:
                i = i.strip("\n").split("\t")
                if i[0] in special:
                    naa.append(i[0])
                    S[0][sample].append(i[1])
                else:
                    ngg.append(i[0])
                    N[0][sample].append(i[1])
            if len(naa) > len(na):
                na = naa
            if len(ngg) > len(ng):
                ng = ngg
            f.close()
    MT = dict()
    out = open(path.replace("/results_htseq-" + mode +"/","/outputs/") + "/htseq-" + mode + "_counts.txt", 'w')
    print >> out, mode + "_id\t" + "\t".join(samples)
    for j in range(len(ng)):
        r = list()
        for sample in samples:
            if not MT.has_key(sample):
                MT[sample] = 0
            if N[0].has_key(sample):
                if j < len(N[0][sample]):
                    r.append(N[0][sample][j])
                    MT[sample] += int(N[0][sample][j])
                else:
                    r.append("NA")
            else:
                r.append("NA")
        print >> out, ng[j] + "\t" + "\t".join(r)
    out.close()
    out = open(path.replace("/results_htseq-" + mode + "/","/outputs/") + "/htseq-" + mode + "_stats.txt", 'w')
    print >> out, "sample_id\t" + "\t".join(na) + "\tTotalCounts\tLink"
    for sample in samples:
        if S[0].has_key(sample) and MT.has_key(sample):
            if len(S[0][sample]) > 0:
                print >> out, sample + "\t" + "\t".join(S[0][sample]) + "\t" + str(MT[sample]) + "\t" +  '<a href="../results_htseq-'+ mode +'/'+sample+'.tab" target="_blank">+</a>'
            else:
                print >> out, sample + "\t" + "\t".join(["0" for k in range(len(na))]) + "\t0\t"
        else:
            print >> out, sample + "\t" + "\t".join(["0" for k in range(len(na))]) + "\t0\t"
    out.close()
    print "  Data found for " + str(len(N[0])) + " of " + str(len(samples)) + " samples"

    print "  Creating annotation file..."
    D, H, L = get_annotation(mode, path.replace("/results_htseq-" + mode + "/", "/config.txt"))
    out = open(path.replace("/results_htseq-" + mode + "/", "/outputs/") + "/htseq-" + mode +"_annotation.txt", 'w')
    print >> out, "\t".join(H)
    k = 0
    for gene in ng:
        if D.has_key(gene):
            print >> out, D[gene]
            k += 1
        else:
            print >> out, "\t".join(["NA" for i in range(len(H))])
    out.close()
    print "  Recovered annotation for " + str(k) + " of " + str(len(ng)) + " genes"
    print "  Computing RPKMs..."
    out = open(path.replace("/results_htseq-" + mode + "/", "/outputs/") + "/htseq-" + mode + "_rpkms.txt", 'w')
    print >> out, "gene_id\t" + "\t".join(samples)
    for j in range(len(ng)):
        r = list()
        for sample in samples:
            if N[0].has_key(sample):
                if len(N[0][sample]) > j:
                    s = str(round(float(N[0][sample][j]) * 1000000000 / (float(L[ng[j]]) * MT[sample]),4))
                else:
                    s = "NA"
                r.append(s)
            else:
                r.append("NA")
        print >> out, ng[j] + "\t" + "\t".join(r)
    out.close()
    return 1