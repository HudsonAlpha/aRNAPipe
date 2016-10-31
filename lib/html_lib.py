import shutil
import os

modules = ["trimgalore", "fastqc", "kallisto", "star", "star-fusion", "picard", "picard_IS", "htseq-gene", "htseq-exon", "varscan", "gatk"]
module_names = {"trimgalore":"",
                "fastqc":"",
                "kallisto":"",
                "star":"",
                "star-fusion":"",
                "picard":"",
                "picard_IS":"",
                "htseq-gene":"",
                "htseq-exon":"",
                "varscan":"",
                "gatk":""}


def get_menu(config, ns):
    enabled = dict()
    for module in modules:
        if config.has_key(module):
            if int(config[module][0].split("/")[0]) > 0:
                enabled[module] = 1
    menu = list()
    menu.append('<h1><a #highlight="" href="./summary.html">Summary</a></h1>')
    menu.append('<h1><a #highlight="" href="./hpc.html">HPC statistics</a></h1>')
    if enabled.has_key("fastqc") or enabled.has_key("trimgalore"):
        menu.append('<h1>Raw QC:</h1>')
    if enabled.has_key("trimgalore"):
        menu.append('<h2><a #highlight="" href="./trim.html">- TrimGalore/Cutadapt</a></h2>')
    if enabled.has_key("fastqc"):
        menu.append('<h2><a #highlight="" href="./fastqc.html">- FastQC</a></h2>')
    if enabled.has_key("star") or enabled.has_key("picard") or enabled.has_key("htseq-gene") or enabled.has_key("htseq-exon"):
        menu.append('<h1>Alignment QC:</h1>')
    if enabled.has_key("picard"):
        menu.append('<h2><a #highlight="" href="./picard.html">- Picard</a></h2>')
    if enabled.has_key("picard_IS"):
        menu.append('<h2><a #highlight="" href="./picard-is.html">- Picard Insert Size</a></h2>')
    if enabled.has_key("star"):
        menu.append('<h2><a #highlight="" href="./star.html">- STAR</a></h2>')
    if enabled.has_key("kallisto"):
        menu.append('<h2><a #highlight="" href="./kallisto.html">- KALLISTO</a></h2>')
    if enabled.has_key("htseq-gene"):
        menu.append('<h2><a #highlight="" href="./htseq-gene.html">- HTseq-Gene</a></h2>')
    if enabled.has_key("htseq-exon"):
        menu.append('<h2><a #highlight="" href="./htseq-exon.html">- HTseq-Exon</a></h2>')
    if ns > 1:
        if enabled.has_key("star") or enabled.has_key("htseq-gene") or enabled.has_key("htseq-exon") or enabled.has_key("kallisto"):
            menu.append('<h1>Count statistics:</h1>')
            menu.append('<h2><a #highlight="" href="./downloads.html">- DOWNLOADS</a></h2>')
        if enabled.has_key("star"):
            menu.append('<h2><a #highlight="" href="./star2.html">- STAR</a></h2>')
        if enabled.has_key("kallisto"):
            menu.append('<h2><a #highlight="" href="./kallisto2.html">- KALLISTO</a></h2>')
        if enabled.has_key("htseq-gene"):
            menu.append('<h2><a  href="./htseq-gene2.html">- HTseq-Gene</a></h2>')
        if enabled.has_key("htseq-exon"):
            menu.append('<h2><a #highlight="" href="./htseq-exon2.html">- HTseq-Exon</a></h2>')
    if enabled.has_key("gatk") or enabled.has_key("varscan"):
        menu.append('<h1>Variant calling:</h1>')
    if enabled.has_key("varscan"):
        menu.append('<h2><a #highlight="" href="./varscan.html">- VARSCAN</a></h2>')
    if enabled.has_key("gatk"):
        menu.append('<h2><a #highlight="" href="./gatk.html">- GATK</a></h2>')
    if enabled.has_key("star-fusion"):
        menu.append('<h1>Gene fusions:</h1>')
        menu.append('<h2><a #highlight="" href="./star-fusion.html">- Star-Fusion</a></h2>')
    menu = "\n".join(menu)
    return menu

def print_samples(path,config):
    analysis = ['trimgalore', 'fastqc', 'kallisto', 'star', 'star-fusion', 'picard', "htseq-gene", "htseq-exon", "picard_IS", "varscan", 'gatk']
    sta= {"trimgalore":"TrimGalore", "fastqc":"FastQC","star":"STAR","star-fusion":"STAR-Fusion","picard":"PicardQC","kallisto":"Kallisto","htseq-gene":"HTseq-gene","htseq-exon":"HTseq-exon", "picard_IS":"Picard-InsertSize", "varscan":"VARSCAN", "gatk":"GATK"}
    # SAMPLES LIST
    samples = dict()
    f   = open(path + "/samples.list",'r')
    hs = f.readline().strip("\n").split("\t")
    for i in f:
        i = i.strip("\n").split("\t")
        if i[0] != "":
            samples[i[0]] = i[1:]
    f.close()
    # scan
    results  = dict()
    for i in analysis:
        if config.has_key(i):
            sok = dict()
            if os.path.exists(path + "/results_" + i + "/samples_ok.txt"):
                f = open(path + "/results_" + i + "/samples_ok.txt", 'r')
                for j in f:
                    sok[j.strip("\n")] = 1
                f.close()
            results[i] = dict()
            if i == "trimgalore":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        filess = y
                        for f in filess:
                            f = f.split("/")[-1]
                            link = "../results_trimgalore/" + f + "_trimming_report.txt"
                            link = '<a href="LINK" target="_blank">OK</a>'.replace("LINK", link)
                            res.append(link)
                    else:
                        res.append("FAIL")
                    results["trimgalore"][x] = " / ".join(res)
            elif i=="fastqc":
                for x, y in sorted(samples.iteritems()):
                    res    = []
                    if sok.has_key(x):
                        filess = y
                        for f in filess:
                            f = f.split("/")[-1]
                            link = "../results_fastqc/"+f.replace(".fastq","").replace(".gz","")+"_fastqc/fastqc_report.html"
                            link = '<a href="LINK" target="_blank">OK</a>'.replace("LINK",link)
                            res.append(link)
                    else:
                        res.append("FAIL")
                    results["fastqc"][x] = " / ".join(res)
            elif i=="star":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        link = "../results_star/" + x + "_Aligned.out.sam"
                        link = '<a href="LINK" target="_blank">BAM-OK</a>'.replace("LINK", link)
                        res.append(link)
                        link = "../results_star/" + x + "_ReadsPerGene.out.tab"
                        link = '<a href="LINK" target="_blank">COUNTS-OK</a>'.replace("LINK", link)
                        res.append(link)
                        link = "../results_star/" + x + "_SJ.out.tab"
                        link = '<a href="LINK" target="_blank">SJ-OK</a>'.replace("LINK", link)
                        res.append(link)
                    else:
                        res.append("BAM-FAIL")
                        res.append("COUNTS-FAIL")
                        res.append("COUNTS-FAIL")
                    results["star"][x] = " / ".join(res)
            elif i=="kallisto":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        link = "../results_kallisto/" + x + "/abundance.tsv"
                        link = '<a href="LINK" target="_blank">OK</a>'.replace("LINK", link)
                        res.append(link)
                    else:
                        res.append("FAIL")
                    results["kallisto"][x] = " / ".join(res)
            elif i=="star-fusion":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        link = "../results_star-fusion/" + x + ".fusion_candidates.txt"
                        link = '<a href="LINK" target="_blank">OK</a>'.replace("LINK", link)
                        res.append(link)
                    else:
                        res.append("FAIL")
                    results["star-fusion"][x] = " / ".join(res)
            elif i=="picard_IS":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        link = "../results_picard_IS/" + x + ".txt"
                        link = '<a href="LINK" target="_blank">OK</a>'.replace("LINK", link)
                        res.append(link)
                    else:
                        res.append("FAIL")
                    results["picard_IS"][x] = " / ".join(res)
            elif i=="picard":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        link = "../results_picard/" + x + ".general.qc"
                        link = '<a href="LINK" target="_blank">GENERAL-OK</a>'.replace("LINK", link)
                        res.append(link)
                        link = "../results_picard/" + x + ".protein_coding.qc"
                        link = '<a href="LINK" target="_blank">PC-OK</a>'.replace("LINK", link)
                        res.append(link)
                        link = "../results_picard/" + x + ".ribosomal.qc"
                        link = '<a href="LINK" target="_blank">RB-OK</a>'.replace("LINK", link)
                        res.append(link)
                    else:
                        res.append("GENERAL-FAIL")
                        res.append("PC-FAIL")
                        res.append("RB-FAIL")
                    results["picard"][x] = " / ".join(res)
            elif i=="htseq-gene":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        link = "../results_htseq-gene/" + x + ".tab"
                        link = '<a href="LINK" target="_blank">OK</a>'.replace("LINK", link)
                        res.append(link)
                    else:
                        res.append("FAIL")
                    results["htseq-gene"][x] = " / ".join(res)
            elif i=="htseq-exon":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        link = "../results_htseq-exon/" + x + ".tab"
                        link = '<a href="LINK" target="_blank">OK</a>'.replace("LINK", link)
                        res.append(link)
                    else:
                        res.append("FAIL")
                    results["htseq-exon"][x] = " / ".join(res)
            elif i=="varscan":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        link = "../results_varscan/" + x + ".vcf"
                        link = '<a href="LINK" target="_blank">VCF</a>'.replace("LINK", link)
                        res.append(link)
                    else:
                        res.append("FAIL")
                    results["varscan"][x] = " / ".join(res)
            elif i=="gatk":
                for x, y in sorted(samples.iteritems()):
                    res = []
                    if sok.has_key(x):
                        link1 = "../results_gatk/" + x + ".vcf"
                        link1 = '<a href="LINK" target="_blank">VCF</a>'.replace("LINK", link1)
                        link2 = "../results_gatk/" + x + ".filt.vcf"
                        link2 = '<a href="LINK" target="_blank">VCF_FILT</a>'.replace("LINK", link2)
                        res.append(link1 + "/" + link2)
                    else:
                        res.append("FAIL")
                    results["gatk"][x] = " / ".join(res)
    n = "<th bgcolor='#A8A8A8'>Sample</th>"
    for i in hs[1:]:
        n = n + "<th bgcolor='#A8A8A8'> Size "+i+" (Gb)</th>"
    for i in range(len(analysis)):
        if results.has_key(analysis[i]):
            n = n +"<th bgcolor='#A8A8A8'>"+sta[analysis[i]]+"</th>"
    thead = "<thead><tr>"+n+"</tr></thead>"
    tab = list()
    for i in sorted(samples.keys()):
        n = ["<td bgcolor='#A8A8A8'>"+i+"</td>"]
        for j in range(len(hs[1:])):
            try:
                n.append("<td bgcolor='#A8A8A8'>" + str(round(os.stat(samples[i][j]).st_size/1000000000.0, 2)) + "</td>")
            except:
                n.append("<td bgcolor='#A8A8A8'>NA</td>")
        for a in analysis:
            if results.has_key(a):
                cl = "#00CC66"
                if "FAIL" in results[a][i]:
                    cl = "#CC3300"
                n.append("<td bgcolor='"+cl+"'>"+results[a][i]+"</td>")
        tab.append("<tr>"+"".join(n)+"</tr>")
    return '<table id="DT" class="display">' + thead + "<tbody>" + "\n".join(tab) + "</tbody></table>"


def config_file(path, fname):
    f   = open(path + "/" + fname,'r')
    n   = list()
    for i in f:
        i = i.strip("\n").split("\t")
        if len(i) >= 2:
            g = "<tr><td bgcolor='#00CC66'>"+i[0]+"</td><td bgcolor='#00CC66'>"+i[1]+"</td></tr>"
        else:
            g = "<tr><td colspan='2' bgcolor='#C0C0C0'>" + i[0] + "</td></tr>"
        n.append(g)
    f.close()
    tab = "<table>"+"".join(n)+"</table>"
    return tab


def print_table_default(datafile, index, select):
    palette = ["#00FA9A", "#AFEEEE", "#D8BFD8", "#DEB887", "#D3D3D3", "#EEE8AA"]
    if not os.path.exists(datafile):
        return ""
    f = open(datafile, 'r')
    h = f.readline().strip("\n").split("\t")
    if len(select) == 0:
        select = range(len(h))
    n = ""
    for i in select:
        n = n + "<th align='center' bgcolor='#A8A8A8'>" + h[i] + "</th>"
    if index < 0:
        n = '<table id="DT" class="display"><thead><tr>' + n + '</tr></thead><tbody>'
    else:
        n = '<table><thead><tr>' + n + '</tr></thead><tbody>'
    M = dict()
    r = 0
    for i in f:
        i = i.strip("\n").split("\t")
        if len(i) > 0:
            temp = ""
            if index > -1:
                if not M.has_key(i[index]):
                    M[i[index]] = palette[r % len(palette)]
                    r += 1
            for k in select:
                j = i[k]
                if (j == "-1") or (j.startswith("NA ") or j.endswith(" NA") or j == "NA"):
                    temp = temp + "<td align='center' bgcolor='#CC3300'>" + j + "</td>"
                else:
                    if index < 0:
                        temp = temp + "<td align='center' bgcolor='#00CC66'>" + j + "</td>"
                    else:
                        temp = temp + "<td align='center' bgcolor='" + M[i[index]] + "'>" + j + "</td>"
            n = n + "<tr>" + temp + "</tr>"
    n += "</tbody></table>"
    return n


def check_project(path):
    print "> Checking project path..."
    if path == "":
        exit("Parameter -p is required.")
    if not os.path.exists(path):
        exit("Path to project folder not found.")
    if path.endswith("/"):
        path = path[0:(len(path)-1)]
    if not os.path.exists(path + "/config.txt"):
        exit("Project configuration file not found: " + path + "/config.txt")
    if not os.path.exists(path + "/samples.list"):
        exit("Project samples file not found: " + path + "/samples.list")
    project = path.split("/")
    project = project[len(project)-1] # project id
    return project, path


def stats_picard(path,samples,config):
    n = os.listdir(path)
    hh = "\t".join(["sample_id","mRNA-Coding", "mRNA-Ribosomal", "mRNA-Others","Intron","Intergenic","Unmapped","Med-CVcov","Med-5bias","Med-3bias","Med-5/3bias"])
    if config.has_key("picard") and ("results_picard" in n):
        files  = os.listdir(path+"/results_picard")
        n      = [".general.qc",".protein_coding.qc",".ribosomal.qc"]
        names  = ["Sample","BASES (N)","ALIGN (N)","ALIGN (%)","INTRON (N)","INTRON (%)",
                  "INTERGEN (N)","INTERGEN (%)","MRNA (N)","MRNA (%)","Protein coding (%)",
                  "Ribosomal (%)","Others (%)","Coding (%)","UTR (%)","MEDIAN CV COVERAGE","MEDIAN 5' BIAS","MEDIAN 3' BIAS","MEDIAN 5'-3' BIAS"]
        g = ["MEDIAN_CV_COVERAGE","MEDIAN_5PRIME_BIAS","MEDIAN_3PRIME_BIAS","MEDIAN_5PRIME_TO_3PRIME_BIAS"]
        out = open(path + "/outputs/stats_picard.txt",'w')
        table  = ['<tr><th align="center" colspan="2" bgcolor="#A8A8A8"></th><th align="center" colspan="2" bgcolor="#A8A8A8">ALIGNED BASES</th><th align="center" colspan="6" bgcolor="#A8A8A8">GENOME LOCATION</th><th align="center" colspan="3" bgcolor="#A8A8A8">MRNA TYPE</th><th align="center" colspan="2" bgcolor="#A8A8A8">MRNA LOCATION</th><th align="center" colspan="4" bgcolor="#A8A8A8">OTHERS</th></tr>']
        header = ""
        for i in names:
            header = header+"<th bgcolor='#A8A8A8'>"+i+"</th>"
        print >> out, hh
        header = "<tr>"+header+"</tr>"
        table.append(header)
        for i in sorted(samples.keys()):
            stats = [{},{},{}]
            heads = list()
            ex = 0
            for k in range(3):
                if i+n[k] in files:
                    f = open(path+"/results_picard"+"/"+i+n[k],'r')
                    nx = 0
                    kdiff0 = 0
                    for ii in f:
                        if ii.startswith("PF_BASES"):
                            head = ii.rstrip().split("\t")
                            nx   = 1
                            for kk in head:
                                if not (kk in heads):
                                    heads.append(kk)
                        elif nx == 1:
                            vals = ii.strip("\n").split("\t")
                            for kk in range(len(head)):
                                stats[k][head[kk]] = vals[kk]
                                if vals[kk] != "":
                                    if float(vals[kk]) > 0:
                                        kdiff0 += 1
                            nx = 0
                    if kdiff0 == 0:
                        ex = 1
                        break
                    f.close()
                else:
                    ex = 1
                    break
            if ex ==1:
                tr = "<td bgcolor='#CC3300'>"+i+"</td>"
                o = i
                for ii in range(18):
                    tr += "<td bgcolor='#CC3300'>NA</td>"
                    o += "\tNA"
                tr = "<tr>"+tr+"</tr>"
                table.append(tr)
                s = ["NA" for ii in range(10)]
            else:
                align  = int(stats[0]["PF_ALIGNED_BASES"])
                cod  = int(stats[0]["CODING_BASES"])
                utr  = int(stats[0]["UTR_BASES"])
                intr = int(stats[0]["INTRONIC_BASES"])
                inter = int(stats[0]["INTERGENIC_BASES"])
                mrna = cod + utr
                pc   = int(stats[1]["CODING_BASES"])+int(stats[1]["UTR_BASES"])
                rb   = int(stats[2]["CODING_BASES"])+int(stats[2]["UTR_BASES"])
                ot   = mrna - pc - rb
                pc   = str(round(100*float(pc)/align,3))
                rb   = str(round(100*float(rb)/align,3))
                ot   = str(round(100*float(ot)/align,3))
                tr = "<td bgcolor='#B8B8B8'>"+i+"</td>"
                tr +="<td bgcolor='#00CC66'>"+stats[0]["PF_BASES"]+"</td>"
                tr +="<td bgcolor='#00CC66'>"+str(align)+"</td>"
                tr +="<td bgcolor='#99CC66'>"+str(round(100*float(align)/float(stats[0]["PF_BASES"]),2))+"</td>"
                tr +="<td bgcolor='#00CC66'>"+str(intr)+"</td>"
                tr +="<td bgcolor='#99CC66'>"+str(round(100*float(intr)/align,2))+"</td>"
                tr +="<td bgcolor='#00CC66'>"+str(inter)+"</td>"
                tr +="<td bgcolor='#99CC66'>"+str(round(100*float(inter)/align,2))+"</td>"
                tr +="<td bgcolor='#00CC66'>"+str(mrna)+"</td>"
                tr +="<td bgcolor='#99CC66'>"+str(round(100*float(mrna)/align,2))+"</td>"
                tr +="<td bgcolor='#99CC66'>"+pc+"</td>"
                tr +="<td bgcolor='#99CC66'>"+rb+"</td>"
                tr +="<td bgcolor='#99CC66'>"+ot+"</td>"
                tr +="<td bgcolor='#99CC66'>"+str(round(100*float(cod)/align,2))+"</td>"
                tr +="<td bgcolor='#99CC66'>"+str(round(100*float(utr)/align,2))+"</td>"
                o = i + "\t" + stats[0]["PF_BASES"] + "\t" + str(align) + "\t" + str(round(100*float(align)/float(stats[0]["PF_BASES"]),2)) + "\t" + str(intr) + "\t" + str(round(100*float(intr)/align,2)) + "\t" + str(inter) + "\t" + str(round(100*float(inter)/align,2)) + "\t" + str(mrna) + "\t" + str(round(100*float(mrna)/align,2)) + "\t" +pc+ "\t" +rb+ "\t" +ot+ "\t" +str(round(100*float(cod)/align,2))+ "\t" + str(round(100*float(utr)/align,2))
                for ix in g:
                    if stats[0][ix]!="?":
                        tr +="<td bgcolor='#99CC66'>"+str(round(float(stats[0][ix]),3))+"</td>"
                        o += "\t" + str(round(float(stats[0][ix]),3))
                    else:
                        tr +="<td bgcolor='#99CC66'>0</td>"
                        o += "\t0"
                tr = "<tr>"+tr+"</tr>"
                table.append(tr)
                o = o.split("\t")
                st= 0
                s = []
                for ind in [10, 11, 12, 5, 7, 3]:
                    s.append(str(round(float(o[ind]) * float(o[2])/float(o[1]),3)))
                    if ind != 3:
                        st += float(o[ind]) * float(o[2])/float(o[1])
                for ind in [15,16,17,18]:
                    if o[ind] != "0":
                        s.append(str(round(float(o[ind]) * 100,3)))
                    else:
                        s.append("0")
                s[5] = str(round(100 - st,3))
            print >> out, i + "\t" + "\t".join(s)
        out.close()
        return "<table>"+"\n".join(table)+"</table>"
    else:
        return ""


def stats_picard_2(path,samples,config):
    n = os.listdir(path)
    hh = "\t".join(['sample_id','MEDIAN_INSERT_SIZE','MEDIAN_ABSOLUTE_DEVIATION','MIN_INSERT_SIZE',
                    'MAX_INSERT_SIZE','MEAN_INSERT_SIZE','STANDARD_DEVIATION','READ_PAIRS', 'LINK_TXT', 'LINK_PDF'])
    print n
    print config
    print hh
    if config.has_key("picard_IS") and ("results_picard_IS" in n):
        files  = os.listdir(path+"/results_picard_IS")
        out = open(path + "/outputs/stats_picard2.txt",'w')
        print >> out, hh
        for i in sorted(samples.keys()):
            if i + ".txt" in files:
                f = open(path+"/results_picard_IS"+"/"+i+".txt",'r')
                k = 0
                while (1):
                    j = f.readline()
                    if j.startswith("MEDIAN_INSERT_SIZE"):
                        j = f.readline().strip("\n").split("\t")
                        break
                    k += 1
                    if k > 10 or len(j) == 0:
                        j = ['NA' for i in range(7)]
                        break
                print >> out, "\t".join([i] + j + ['<a href="../results_picard_IS/' + i + '.txt" target="_blank">+</a>', '<a href="../results_picard_IS/' + i + '.pdf" target="_blank">+</a>'])
                f.close()
            else:
                print >> out, "\t".join([i] + ['NA' for i in range(9)])
        out.close()
    return 1

def skeleton(path, path2html):
    print "> Building HTML and OUTPUT folders skeletons..."
    print "  - Path: " + path
    print "  - Libs: " + path2html
    # Creates output directiories
    if os.path.exists(path + "/outputs"):
        os.system("rm -r " + path + "/outputs")
    os.mkdir(path + "/outputs")
    # Creates HTML directories
    n = os.listdir(path)
    if "HTML" in n:
        os.system("rm -r " + path + "/HTML")
    if not ("HTML" in os.listdir(path)):
        os.mkdir(path + "/HTML")
    os.mkdir(path + "/HTML/html")
    # Copy lightbox and jquery
    shutil.copy(path2html + "/html/style.css", path + "/HTML/html/style.css")
    shutil.copy(path2html + "/html/lytebox.js", path + "/HTML/html/lytebox.js")
    shutil.copy(path2html + "/html/lytebox.css", path + "/HTML/html/lytebox.css")
    shutil.copy(path2html + "/html/jquery-1.12.0.js", path + "/HTML/html/jquery-1.12.0.js")
    shutil.copy(path2html + "/html/jquery.dataTables.min.js", path + "/HTML/html/jquery.dataTables.min.js")
    shutil.copy(path2html + "/html/jquery.dataTables.min.css", path + "/HTML/html/jquery.dataTables.min.css")
    shutil.copy(path2html + "/html/dataTables.colReorder.min.js", path + "/HTML/html/dataTables.colReorder.min.js")
    shutil.copy(path2html + "/html/amcharts.js", path + "/HTML/html/amcharts.js")
    shutil.copy(path2html + "/html/serial.js", path + "/HTML/html/serial.js")
    shutil.copy(path2html + "/html/xy.js", path + "/HTML/html/xy.js")
    os.system("cp -r " + path2html + "/html/images " + path + "/HTML/html/")


# def check_samples(path):
#     # Parses the samples file
#     print "> Parsing samples file..."
#     try:
#         f = open(path + "/samples.list", 'r')
#         i = f.readline().strip("\n").split("\t")
#         index = [-1,-1,-1]
#         if ("FASTQ_1" in i) and ("FASTQ_2" in i):
#             k = "paired-end"
#             for r in range(len(i)):
#                 if i[r]=="SampleID":
#                     index[0] = r
#                 if i[r]=="FASTQ_1":
#                     index[1] = r
#                 if i[r]=="FASTQ_2":
#                     index[2] = r
#         else:
#             k = "single-end"
#             for r in range(len(i)):
#                 if i[r]=="SampleID":
#                     index[0] = r
#                 if i[r]=="FASTQ":
#                     index[1] = r
#         samples = dict()
#         for i in f:
#             i = i.strip("\n").split("\t")
#             if len(i)>1:
#                 if k=="paired-end":
#                     samples[i[index[0]]] = [[i[index[1]],i[index[2]]],["Type",k]]
#                 else:
#                     samples[i[index[0]]] = [[i[index[1]]],["Type",k]]
#         f.close()
#         return samples
#     except:
#         exit("Error checking samples file: " + path + "/samples.list")


def build_amcharts(input, output, prog, pname, path, html_table, project, lmenu):
    out = open(output, 'w')
    f   = open(input, 'r')
    it = ""
    for i in f:
        i = i.strip("\n")
        if i.startswith("#ITERATOR"):
            pattern = i.replace("#ITERATOR=","")
            f2 = open(path + "/outputs/" + pname[1] + "_pca.txt", 'r')
            h  = f2.readline()
            it = list()
            for j in f2:
                np = pattern
                j = j.strip("\n").split("\t")
                for k in range(len(j)):
                    np = np.replace("#VAR"+str(len(j)-k-1),j[len(j)-k-1])
                it.append(np)
            it = ", ".join(it)
            f2.close()
        else:
            if prog + "2.html" in i:
                i = i.replace("#HIGHLIGHT", 'style="color:#808080"')
            print >> out, i.replace("#LATMENU",lmenu).replace("#PROG", pname[0]).replace("#PROJECT", project).replace("#SITERATOR", it).replace("#HIGHTLIGHT", "").replace("#TABLE", html_table)
    f.close()
    out.close()


def check_config(path):
    # Parses the configuration file
    print "> Parsing configuration file..."
    try:
        z = ["trimgalore", "fastqc", "star", "star-fusion", "picard", "htseq-gene", "htseq-exon", "kallisto", "picard_IS", "varscan", "gatk"]
        f = open(path + "/config.txt", 'r')
        analysis = dict()
        analysis["cluster"] = dict()
        analysis["programs"] = dict()
        for i in f:
            if not i.startswith("#"):
                i = i.strip("\n").split("\t")
                if len(i) > 1:
                    if i[0] in z:
                        if int(i[1].split("/")[0]) > 0:
                            analysis[i[0]] = [i[1], "results_" + i[0], dict()]
                    elif i[0] in ["wt", "q"]:
                        analysis["cluster"][i[0]] = i[1]
                    elif i[0] == "star_args_own":
                        if analysis["programs"]["star_args"] == "own":
                            analysis["programs"]["star_args"] = i[1]
                    elif i[0] == "starfusion_own":
                        if analysis["programs"]["starfusion"] == "own":
                            analysis["programs"]["starfusion"] = i[1]
                    elif i[0] == "genome_build":
                        analysis[i[0]] = i[1]
                    else:
                        analysis["programs"][i[0]] = i[1]
        f.close()
        return analysis
    except:
        exit("Error checking configuration file: " + path + "/config.txt")


# def print_config(config,path):
#     table = list()
#     table.append(["Analysis", "Processors", "Folder", "Timestamp",
#                   "TStart","TEnd","Success","CPU-Time","MaxMemo",
#                   "AveMemo","MaxSwap","Parameters"])
#     for i in modules:
#         if config.has_key(i):
#             n = check_log_cluster(path, config[i][1])
#             st = []
#             if len(config[i][2]) > 0:
#                 for v,w in config[i][2].iteritems():
#                     st.append(v+": "+w)
#             st = "<br>".join(st)
#             for j in range(len(n)):
#                 tt = [module_names[i]]+[config[i][0],"./"+config[i][1]]+n[j]+[st]
#                 table.append(tt)
#     n = ""
#     for i in table[0]:
#         n = n+"<th bgcolor='#A8A8A8'>"+i+"</th>"
#     n = ["<tr>"+n+"</tr>"]
#     for i in table[1:]:
#         temp = ""
#         for j in i:
#             if "NA" in j:
#                 temp = temp+"<td bgcolor='#CC3300'>"+j+"</td>"
#             else:
#                 temp = temp+"<td bgcolor='#00CC66'>"+j+"</td>"
#         n.append("<tr>"+temp+"</tr>")
#     return n


# def check_log_cluster(path, val):
#     t = os.listdir(path)
#     if not (val in t):
#         return ["NA","NA","NA","NA","NA","NA","NA"]
#     t  = os.listdir(path + "/" + val)
#     t2 = list()
#     for i in t:
#         if i.startswith("log_cluster_") and (("scheduler" in i) == False):
#             if len(i.split("_")) >= 4 :
#                 t2.append(i)
#     if len(t2) == 0:
#         return [["NA","NA","NA","NA","NA","NA","NA","NA"]]
#     n = list()
#     for jv in sorted(t2):
#         f = open(path + "/" + val + "/" + jv,'r')
#         ts = ""
#         te = ""
#         suc= "No"
#         cpu_time = ""
#         max_memo = ""
#         ave_memo = ""
#         max_swap = ""
#         pid      = "_".join(jv.split("_")[2:4]).replace(".txt","")
#         for i in f:
#             if i.startswith("Started at"):
#                 ts = i.rstrip().replace("Started at ","")
#             if i.startswith("Results reported on"):
#                 te = i.rstrip().replace("Results reported on ","")
#             if i.startswith("Successfully completed."):
#                 suc = "Yes"
#             if "CPU time :" in i:
#                 cpu_time = " ".join(i.rstrip().split()[3:5])
#             if "Max Memory :" in i:
#                 max_memo = " ".join(i.rstrip().split()[3:5])
#             if "Average Memory :" in i:
#                 ave_memo = " ".join(i.rstrip().split()[3:5])
#             if "Max Swap :" in i:
#                 max_swap = " ".join(i.rstrip().split()[3:5])
#         f.close()
#         n.append([pid,ts,te,suc,cpu_time,max_memo,ave_memo,max_swap])
#     return n


def bar_getdata (filename, head, cols_bar, cols_line):
    # LOAD DATA
    if not os.path.exists(filename):
        return ""
    f = open(filename, 'r')
    h = f.readline().strip("\n").split("\t")
    if (len(cols_bar)==0) and (len(cols_line)==0):
        cols_bar = range(1, len(h))
    D = list()
    for i in f:
        i = i.strip("\n").split("\t")
        n = list()
        for j in range(len(h)):
            if j == head:
                n.append('"' + h[j] + '": "' + i[j] + '"')
            elif (j in cols_bar) or (j in cols_line):
                if i[j] != "NA":
                    n.append('"' + h[j] + '": ' + i[j])
                else:
                    n.append('"' + h[j] + '": 0')
        D.append("{" + ", ".join(n) + "}")
    D = "var chartData = [" + ", ".join(D) + "];"
    return D


def build_from_template(prog, project, data, html_table, html_table2, output, template, lmenu):
    out = open(output,'w')
    print output
    f = open(template, 'r')
    print template
    r = output.split("/")[-1]
    print r
    for i in f:
        i = i.strip("\n")
        print i
        if r in i:
            i = i.replace("#HIGHLIGHT", 'style="color:#808080"')
        print >> out, i.replace("#LATMENU",lmenu).replace("#PROG", prog).replace("#PROJECT", project).replace("#DATA", data).replace("#HIGHTLIGHT", "").replace("#TABLE2", html_table2).replace("#TABLE", html_table)
    f.close()
    out.close()