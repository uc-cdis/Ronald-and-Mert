import HTSeq
import re
import os

baseDirPath = "/mnt/reference_bams/"
os.chdir(baseDirPath)
bamFileRegex = re.compile(r'GRCh37\.HumanBodyMap\.(.*)\.1\.bam')
gtf_file = HTSeq.GFF_Reader("Homo_sapiens.GRCh37.70.gtf", end_included = True)
for root, dirs, files in os.walk('.'):
    for f in files:
        regexSearch = bamFileRegex.search(f)
        if regexSearch != None:
            tissueSite = regexSearch.group(1)
        if not os.path.isfile(tissueSite + ".counts"):
            print f
            outfile = open(tissueSite + ".counts", 'w')
            bam_reader = HTSeq.BAM_Reader(f)
            exons = HTSeq.GenomicArrayOfSets("auto", stranded = False)
            for feature in gtf_file:
                if feature.type == "exon":
                    exons[feature.iv] += feature.name
            counts = {}
            for feature in gtf_file:
                if feature.type == "exon":
                    counts[feature.name] = 0
            for alnmt in bam_reader:
                if alnmt.aligned:
                    iset = None
                    for iv2, step_set in exons[alnmt.iv].steps():
                        if iset is None:
                            iset = step_set.copy()
                        else:
                            iset.intersection_update(step_set)
                    if len(iset) == 1:
                        counts[list(iset)[0]] += 1
            for name in sorted(counts.keys()):
                outfile.write(name + "\t" + str(counts[name]) + "\n")
            outfile.close()
