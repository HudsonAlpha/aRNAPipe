import os


def fasta2dict(fasta_file, samtools, picard):
    try:
        print "  - genome.dict"
        os.system("java -jar " + picard + "/CreateSequenceDictionary.jar R=" + fasta_file + " O=" + fasta_file.replace(".fa", ".dict") + " &>/dev/null")
        print "  - genome.fa.fai"
        os.system(samtools + " faidx " + fasta_file)
        return 1
    except:
        return 0

def refflat(path_gtf2gp, output_dir):
    try:
        print "  - Global file: 'genesets.refFlat'"
        outa = open(output_dir + "/script.sh", 'w')
        temp = output_dir + "/temp.rf"
        command1 = path_gtf2gp + " -genePredExt -geneNameAsName2 " + output_dir + "/genesets.gtf" + " " + temp
        command2 = "paste <(cut -f 12 " + temp +") <(cut -f 1-10 " + temp +") > " + temp.replace("temp.rf", "genesets.refFlat")
        print >> outa, command1
        print >> outa, command2
        print >> outa, "rm " + temp
        ofold = output_dir + "/refFlats"
        types = dict()
        if not os.path.exists(ofold):
            os.mkdir(ofold)
        f = open(output_dir + "/genesets.gtf", 'r')
        i = ""
        while 1:
            i = f.readline()
            if not i.startswith("#"):
                break
        if len(i) > 0:
            i = i.strip("\n").split("\t")
            t = i[8].split("gene_biotype ")[1].split(";")[0].replace('"','')
            if not types.has_key(t):
                types[t] = 0
            types[t] += 1
            for i in f:
                i = i.strip("\n").split("\t")
                t = i[8].split("gene_biotype ")[1].split(";")[0].replace('"','')
                if not types.has_key(t):
                    types[t]=0
                types[t]+=1
        f.close()
        N = 0
        for i,j in types.iteritems():
            N += j
        out = dict()
        for i,j in sorted(types.iteritems()):
            print "  - Gene\t"+i+"\t"+str(j)+"\t"+str(round(100*float(j)/N,2))+"%"
            if j > 0:
                out[i] = open(ofold+"/"+i+".gtf",'w')

        f = open(output_dir + "/genesets.gtf",'r')
        for i in f:
            if i.startswith("#"):
                for j in out.keys():
                    print >> out[j], i.strip("\n")
            else:
                i = i.strip("\n")
                t = i.split("\t")[8].split("gene_biotype ")[1].split(";")[0].replace('"','')
                if out.has_key(t):
                    print >> out[t], i
        f.close()
        for i in out.keys():
            out[i].close()
        n    = os.listdir(ofold)
        temp = output_dir + "/temp.rf"
        for i in n:
            if i.endswith(".gtf"):
                command1 = path_gtf2gp + " -genePredExt -geneNameAsName2 " + ofold + "/" + i + " " + temp
                command2 = "paste <(cut -f 12 " + temp +") <(cut -f 1-10 " + temp +") > " + ofold + "/" + i.replace(".gtf",".refFlat")
                print >> outa, command1
                print >> outa, command2
                print >> outa, "rm " + temp
                print "  - " + i.replace(".gtf",".refFlat")
        outa.close()
        os.system("bash " + output_dir + "/script.sh")
        os.system("rm " + output_dir + "/script.sh")
        return 1
    except:
        return 0

def annotate_gtf(filename):
    try:
        features    = {"gene", "transcript", "exon"}
        transc2gene = dict() # transcript id to gene id for annotation
        exon2gene   = dict() # exon id to gene id for annotation
        loci        = dict() # chrom, start, end
        data = dict()
        for i in features:
            data[i] = dict()
        f = open(filename, 'r')
        for i in f:
            if not i.startswith("#"):
                i = i.strip("\n").split("\t")
                if len(i) == 9:
                    if i[2] in features: # gene/transcript/exon
                        j = i[8]
                        gid = j.split('gene_id "')[1].split('"')[0]  # associated gene id
                        if i[2]+'_id "' in j:
                            feat_id = j.split(i[2]+'_id "')[1].split('"')[0] # feature id
                        else:
                            if i[2]=='exon':
                                feat_id = gid + '_' + j.split('transcript_id "')[1].split('"')[0] + '_' + j.split('exon_number "')[1].split('"')[0]
                        loci[feat_id] = [i[0], i[3], i[4]] # chrom, start, end
                        L = [float(i[3]), float(i[4])]
                        data[i[2]][feat_id] = list()
                        if i[2] == "exon":
                            tid = j.split('transcript_id "')[1].split('"')[0]
                            if not data["gene"].has_key(gid):
                                data['gene'][gid] = list()
                            data["gene"][gid].append(L)
                            if not data["transcript"].has_key(tid):
                                data['transcript'][tid] = list()
                            data["transcript"][tid].append(L)
                            data["exon"][feat_id].append(L)
                            exon2gene[feat_id] = gid
                        elif i[2] == "transcript":
                            transc2gene[feat_id] = gid
        f.close()
        results = dict()
        for i in features:
            results[i] = dict()
        for feat_type, feat_data in data.iteritems():
            for feat_id, exons in feat_data.iteritems():
                ints = list([exons[0]])
                if len(exons) > 1:
                    # Add exons to intervals
                    for exon in exons[1:]:
                        kint = 0
                        for i in range(len(ints)):
                            if (exon[0] <= ints[i][1]) and (ints[i][0] <= exon[1]):
                                # intersection
                                if exon[0] < ints[i][0]:
                                    ints[i][0] = exon[0]
                                if exon[1] > ints[i][1]:
                                    ints[i][1] = exon[1]
                                kint = 1
                                break
                        if kint == 0:
                            ints.append(exon)
                nrem = -1
                while 1:
                    if nrem >= 0:
                        del ints[nrem]
                    nrem = -1
                    if len(ints) > 1:
                        for i in range(len(ints)-1):
                            for j in range(i+1, len(ints)):
                                if (ints[i][0] <= ints[j][1]) and (ints[j][0] <= ints[i][1]):
                                    # intersection
                                    if ints[j][0] < ints[i][0]:
                                        ints[i][0] = ints[j][0]
                                    if ints[j][1] > ints[i][1]:
                                        ints[i][1] = ints[j][1]
                                    nrem = j
                    else:
                        break
                    if nrem == -1:
                        break
                n = 0
                for i in ints:
                    n += i[1]-i[0]+1
                results[feat_type][feat_id] = str(n)
        for feat_type, feat_data in results.iteritems():
            print "  - " + filename.replace(".gtf", "." + feat_type + ".txt")
            out = open(filename.replace(".gtf", "." + feat_type + ".txt"), 'w')
            print >> out, "gene\ttranscript\texon\tchrom\tstart\tend\tlength"
            for feat_id, feat_length in sorted(feat_data.iteritems()):
                if loci.has_key(feat_id):
                    if feat_type == "gene":
                        print >> out, feat_id + "\t.\t.\t" + "\t".join(loci[feat_id]) + "\t" + feat_length
                    elif feat_type == "transcript":
                        print >> out, transc2gene[feat_id] + "\t" + feat_id + "\t.\t" +  "\t".join(loci[feat_id]) + "\t" + feat_length
                    elif feat_type == "exon":
                        print >> out, exon2gene[feat_id] + "\t.\t" + feat_id + "\t" +  "\t".join(loci[feat_id]) + "\t" + feat_length
            out.close()
        return 1
    except:
        return 0






