import os
import sys
import optparse
import lib.vcrparser as vcrparser
import lib.spider_stats as spider_stats
import lib.html_lib as html

#########################################################################
# PARSER
#########################################################################
desc = "aRNApipe: SPIDER module"
parser = optparse.OptionParser(description = desc)
parser.add_option("-p", "--path", dest = "path", default = "", help = "Required: Path to the project folder")
(opt, args) = parser.parse_args()

#########################################################################
# INITIAL ARRANGEMENTS
#########################################################################
pathscript = os.path.dirname(sys.argv[0]) + "/R/" # PATH TO R SCRIPTS
project, path = html.check_project(opt.path) # CHECKS IF PROJECT FOLDER, SAMPLES FILE AND CONFIGURATION FILE EXIST
html.skeleton(path, os.path.dirname(sys.argv[0])) # CREATES THE SKELETON FOR THE HTML AND OUTPUT RESULTS
config = html.check_config(path) # PARSES THE CONFIGURATION FILE
try:
    samples = vcrparser.get_samples(path.replace(project, ""), project, path + "/samples.list") # PARSES THE SAMPLE FILE AND ASSOCIATED FASTQ FILES
except:
    samples = vcrparser.get_samples(path.replace(project, ""), project, path + "/samples.list", no_check=True)
f = open(path + "/samples.list", 'r')
h = f.readline()
samples_ordered = list()
for i in f:
    if len(i.split('\t')) > 0:
        samples_ordered.append(i.split("\t")[0])
f.close()
lmenu = html.get_menu(config, len(samples))

#########################################################################
# PROCESSES AND ARRANGES LOGS AND OUTPUT DATA FROM STAR, KALLISTO AND HTSEQ
#########################################################################
try:
    spider_stats.stats_varcall(path, "varscan")
except:
    print "Unexpected error generating stats of VARSCAN."
try:
    spider_stats.stats_varcall(path, "gatk")
except:
    print "Unexpected error generating stats of GATK."
try:
    spider_stats.stats_trimgalore(path)
except:
    print "Unexpected error generating stats of TRIMGALORE."
try:
    spider_stats.stats_fusion(path)
except:
    print "Unexpected error generating stats of Star-Fusion."
try:
    spider_stats.stats_log(path + "/logs/") # PARSES THE ASSOCIATED LSF LOG FILES AND WRTIES THE RELEVANT PARAMETERS IN OUTPUTS
except:
    print "Unexpected error generating stats of HPC."
try:
    spider_stats.stats_kallisto(path + "/results_kallisto/", samples_ordered) # GENERATES STATISTICS, ESTIMATED COUNTS, TPM AND ANNOTATION MATRICES
except:
    print "Error: Statistics of KALLISTO not yet ready or unexpected error."
try:
    spider_stats.stats_star(path + "/results_star/", samples_ordered)  # GENERATES STATISTICS, COUNTS, RPKM AND ANNOTATION MATRICES
except:
    print "Error: Statistics of STAR not yet ready or unexpected error."
try:
    spider_stats.stats_htseq(path + "/results_htseq-gene/", samples_ordered, "gene") # GENERATES STATISTICS, COUNTS, RPKM AND ANNOTATION MATRICES
except:
    print "Error: Statistics of HTSEQ-GENE not yet ready or unexpected error."
try:
    spider_stats.stats_htseq(path + "/results_htseq-exon/", samples_ordered, "exon") # GENERATES STATISTICS, COUNTS, RPKM AND ANNOTATION MATRICES
except:
    print "Error: Statistics of HTSEQ-EXON not yet ready or unexpected error."
#########################################################################
# HTML SUMMARY WEBPAGE
#########################################################################
try:
    print "> Generating webpage with samples list and configuration..."
    print "  - " + path + "/HTML/summary.html"
    html_table = html.print_samples(path,config) # PROVIDES HTML TABLE WITH SAMPLES STATS
    html_table2 = html.config_file(path, "config.txt") # PROVIDES HTML TABLE WITH CONFIGURATION SETTINGS
    html.build_from_template("SUMMARY", project, "", html_table, html_table2, path+"/HTML/summary.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_SUMMARY.html", lmenu)
except:
    print "  - Not ready"

#########################################################################
# HTML DOWNLOADS
#########################################################################
try:
    print "> Generating webpage with download links..."
    print "  - " + path + "/HTML/downloads.html"
    html.build_from_template("DOWNLOADS", project, "", "", "", path+"/HTML/downloads.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_DOWNLOAD.html", lmenu)
except:
    print "  - Not ready"

#########################################################################
# HPC STATS
#########################################################################
try:
    print "> Generating webpage with HPC and LOG statistics..."
    print "  - " + path + "/HTML/hpc.html"
    os.system("Rscript "+pathscript+"/log_stats.R " + path + "/outputs/log_stats.txt") # PLOT OF HPC USAGE
    html_table = html.print_table_default(path + "/outputs/log_stats.txt", 2, []) # PROVIDES HTML TABLE WITH HPC STATS
    fils = os.listdir(path + "/outputs/")
    gg = ""
    for i in fils:
        if i.endswith(".png") and i.startswith("log_stats_"):
            gg = gg + '<tr bgcolor="#A8A8A8"><td><center><b>Analysis run number '+i.split("_")[2].replace(".png","")+'</b></center></td></tr><tr bgcolor="#00CC66"><td><img src="../outputs/'+i+'" style="width:800px;"></td></tr>'
    html_table2 = "<table>" + gg + "</table>"
    html.build_from_template("HPC", project, "", html_table, html_table2, path+"/HTML/hpc.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_HPC.html", lmenu)
except:
    print "  - Not ready"

#########################################################################
# TRIM_GALORE
#########################################################################
#try:
if config.has_key("trimgalore"):
    print "> Generating webpage with TrimGalore/Cutadapt statistics..."
    print "  - " + path + "/HTML/trim.html"
    html_table = html.print_table_default(path + "/outputs/stats_trim.txt", -1, []) # PROVIDES HTML TABLE WITH HPC STATS
    data = html.bar_getdata(path + "/outputs/stats_trim_plot.txt",0,[],[])
    html.build_from_template("TrimGalore", project, data, html_table, "", path+"/HTML/trim.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_TRIMG.html", lmenu)
# except:
#     print "  - Not ready"

#########################################################################
# FASTQ
#########################################################################
try:
    if config.has_key("fastqc"):
        print "> Generating webpage with fastqc statistics..."
        print "  - " + path + "/HTML/fastqc.html"
        html_table = spider_stats.stats_fastq(path,samples,config) # PROVIDES HTML TABLE WITH FASTQ STATS
        html.build_from_template("FASTQC", project, "", html_table, "", path+"/HTML/fastqc.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_FASTQC.html", lmenu)
except:
    print "  - Not ready"

#########################################################################
# PICARD
#########################################################################
try:
    if config.has_key("picard"):
        print "> Generating webpage with picard statistics..."
        print "  - " + path + "/HTML/picard.html"
        html_table = html.stats_picard(path,samples,config) # PROVIDES HTML TABLE WITH PICARD STATS
        data = html.bar_getdata (path + "/outputs/stats_picard.txt",0,range(1,7), range(7,11))
        html.build_from_template("PICARD", project, data, html_table, "", path+"/HTML/picard.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_PICARD.html", lmenu)
    if config.has_key("picard_IS"):
        print "> Generating webpage with picard insert size statistics..."
        print "  - " + path + "/HTML/picard-is.html"
        x = html.stats_picard_2(path,samples,config)
        html_table = html.print_table_default(path + "/outputs/stats_picard2.txt", -1, [])
        data = html.bar_getdata (path + "/outputs/stats_picard2.txt",0,range(1,2),[])
        html.build_from_template("PICARD-InsertSize", project, data, html_table, "", path+"/HTML/picard-is.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_PICARDIS.html", lmenu)
except:
    print "  - Not ready"

#########################################################################
# HTML STAR-FUSION WEBPAGE
#########################################################################
try:
    if config.has_key("star-fusion"):
        print "> Generating webpage with Star-Fusion results..."
        print "  - " + path + "/HTML/star-fusion.html"
        html_table = html.print_table_default(path + "/outputs/starfusion_aggregate.txt", -1, []) # PROVIDES HTML TABLE WITH HPC STATS
        html.build_from_template("STAR-FUSION", project, "", html_table, "", path+"/HTML/star-fusion.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_STARFUSION.html", lmenu)
except:
    print "  - Not ready"

#########################################################################
# KALLISTO_QC
#########################################################################
try:
    if config.has_key("kallisto"):
        print "> Generating webpage with Kallisto statistics..."
        print "  - " + path + "/HTML/kallisto.html"
        html_table = html.print_table_default(path + "/outputs/kallisto_stats_est_counts.txt", -1, []) # PROVIDES HTML TABLE WITH HPC STATS
        data = html.bar_getdata (path + "/outputs/kallisto_stats_est_counts.txt",0,range(1,2),[])
        html.build_from_template("KALLISTO", project, data, html_table, "", path+"/HTML/kallisto.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_KALLISTO.html", lmenu)
except:
    print "  - Not ready"

#########################################################################
# STAR_QC
#########################################################################
try:
    if config.has_key("star"):
        print "> Generating webpage with STAR statistics..."
        print "  - " + path + "/HTML/star.html"
        html_table = html.print_table_default(path + "/outputs/star_unstranded_stats.txt", -1, []) # PROVIDES HTML TABLE WITH HPC STATS
        data = html.bar_getdata (path + "/outputs/star_unstranded_stats.txt",0,range(1,6),[])
        html.build_from_template("STAR", project, data, html_table, "", path+"/HTML/star.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_STAR.html", lmenu)
except:
    print "  - Not ready"

#########################################################################
# HTSEQ_QC
#########################################################################
try:
    for ij in ["htseq-gene", "htseq-exon"]:
        if config.has_key(ij):
            print "> Generating webpage with "+ij+" statistics..."
            print "  - " + path + "/HTML/"+ij+".html"
            html_table = html.print_table_default(path + "/outputs/"+ij+"_stats.txt", -1, []) # PROVIDES HTML TABLE WITH HPC STATS
            data = html.bar_getdata (path + "/outputs/"+ij+"_stats.txt",0,range(1,7),[])
            html.build_from_template(ij.upper(), project, data, html_table, "", path+"/HTML/" + ij + ".html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_HTSEQ.html", lmenu)
except:
    print "  - Not ready"

#########################################################################
# STAR & HTSEQ STATS ON COUNTS/RPKMS
#########################################################################
try:
    if len(samples) > 1:
        if config["programs"]["strandedness"] == "yes":
            n = {"star":["STAR","star_stranded"],"htseq-gene":["HTseq-count Gene", "htseq-gene"],"htseq-exon":["HTseq-count Exon", "htseq-exon"], 'kallisto': ['Kallisto', 'kallisto']}
        elif config["programs"]["strandedness"] == "no":
            n = {"star":["STAR","star_unstranded"],"htseq-gene":["HTseq-count Gene", "htseq-gene"],"htseq-exon":["HTseq-count Exon", "htseq-exon"], 'kallisto': ['Kallisto', 'kallisto']}
        elif config["programs"]["strandedness"] == "reverse":
            n = {"star":["STAR","star_stranded-reverse"],"htseq-gene":["HTseq-count Gene", "htseq-gene"],"htseq-exon":["HTseq-count Exon", "htseq-exon"], 'kallisto': ['Kallisto', 'kallisto']}
        for prog, pname in n.iteritems():
            if config.has_key(prog):
                os.system("Rscript "+pathscript+"/stats_algs.R " + path + "/outputs/ " + pname[1]) # PLOT OF HPC USAGE
                html_table = html.print_table_default(path + "/outputs/" + pname[1] + "_pca.txt", -1, [0, 1, 2, 3, 4, 6, 7, 8, 12, 13, 14, 15, 17, 18, 19])
                if prog != 'kallisto':
                    html.build_amcharts(os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_PCA.html", path + "/HTML/" + prog + "2.html", prog, pname, path, html_table, project, lmenu)
                else:
                    html.build_amcharts(os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_PCA2.html", path + "/HTML/" + prog + "2.html", prog, pname, path, html_table, project, lmenu)
except:
    print "  - Not ready"

#########################################################################
# VARSCAN
#########################################################################
try:
    if config.has_key("varscan"):
        print "> Generating webpage with VARSCAN statistics..."
        print "  - " + path + "/HTML/varscan.html"
        html_table = html.print_table_default(path + "/outputs/stats_varscan.txt", -1, []) # PROVIDES HTML TABLE WITH HPC STATS
        data = html.bar_getdata (path + "/outputs/stats_varscan.txt",0,[],[])
        html.build_from_template("VARSCAN", project, data, html_table, "", path+"/HTML/varscan.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_GATK.html", lmenu)
except:
    print "  - Not ready"

#########################################################################
# GATK
#########################################################################
try:
    if config.has_key("gatk"):
        print "> Generating webpage with GATK statistics..."
        print "  - " + path + "/HTML/gatk.html"
        html_table = html.print_table_default(path + "/outputs/stats_gatk.txt", -1, []) # PROVIDES HTML TABLE WITH HPC STATS
        data = html.bar_getdata (path + "/outputs/stats_gatk.txt",0,[],[])
        html.build_from_template("GATK", project, data, html_table, "", path+"/HTML/gatk.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_GATK.html", lmenu)
except:
    print "  - Not ready"

#########################################################################
# jSplice
#########################################################################
try:
    if config.has_key("jsplice"):
        print "> Generating webpage with picard statistics..."
        print "  - " + path + "/HTML/jsplice.html"
        html_table = html.print_table_default(path + "/results_jsplice/jSplice_results.txt", -1, []) # PROVIDES HTML TABLE WITH HPC STATS
        html.build_from_template("jSplice", project, "", html_table, "", path+"/HTML/jsplice.html", os.path.dirname(sys.argv[0]) + "/template/TEMPLATE_JSPLICE.html", lmenu)
except:
    print "  - Not ready"
