import os


headers = ["Sample", "MLST", "abcZ", "adk", "aroE", "fumC", "gdh", "pdhC", "pgm", "NgSTAR", "penA NgSTAR", "penA comment",
           "mtrR NgSTAR", "mtrR comment", "porB NgSTAR", "porB comment", "ponA NgSTAR", "ponA comment",
           "gyrA NgSTAR", "gyrA comment", "parC NgSTAR", "parC comment", "23S NgSTAR", "23S comment", "NgMAST",
           "porB NgMAST", "tbpB", "rplF", "rplF species", "rplf species comment", "ppnG coverage", "ppnG depth"]

outstring = snakemake.params.sample
with open(snakemake.input.mlst) as f:
    contig, scheme, profile, abcZ, adk, aroE, fumC, gdh, pdhC, pgm = f.readline().rstrip().split("\t")
outstring += "\t" + profile
outstring += "\t" + abcZ.split('(')[1].split(')')[0]
outstring += "\t" + adk.split('(')[1].split(')')[0]
outstring += "\t" + aroE.split('(')[1].split(')')[0]
outstring += "\t" + fumC.split('(')[1].split(')')[0]
outstring += "\t" + gdh.split('(')[1].split(')')[0]
outstring += "\t" + pdhC.split('(')[1].split(')')[0]
outstring += "\t" + pgm.split('(')[1].split(')')[0]

with open(snakemake.input.ngstar) as f:
    contig, scheme, profile, penA, mtrR, porB, ponA, gyrA, parC, rna23S = f.readline().rstrip().split("\t")
comment_dict = {"penA":{"-":"-"}, "mtrR":{"-":"-"}, "porB":{"-":"-"}, "ponA":{"-":"-"}, "gyrA":{"-":"-"}, "parC":{"-":"-"}, "rna23S":{"-":"-"}}
with open(os.path.join(snakemake.params.mlst_dir, "db", "pubmlst", "ngstar", "NEIS1753.tfa.comments")) as f:
    for line in f:
        allele, comment = line.rstrip().split("\t")
        comment_dict["penA"][allele] = comment
with open(os.path.join(snakemake.params.mlst_dir, "db", "pubmlst", "ngstar", "mtrR.tfa.comments")) as f:
    for line in f:
        allele, comment = line.rstrip().split("\t")
        comment_dict["mtrR"][allele] = comment
with open(os.path.join(snakemake.params.mlst_dir, "db", "pubmlst", "ngstar", "NG_porB.tfa.comments")) as f:
    for line in f:
        allele, comment = line.rstrip().split("\t")
        comment_dict["porB"][allele] = comment
with open(os.path.join(snakemake.params.mlst_dir, "db", "pubmlst", "ngstar", "NG_ponA.tfa.comments")) as f:
    for line in f:
        allele, comment = line.rstrip().split("\t")
        comment_dict["ponA"][allele] = comment
with open(os.path.join(snakemake.params.mlst_dir, "db", "pubmlst", "ngstar", "NG_gyrA.tfa.comments")) as f:
    for line in f:
        allele, comment = line.rstrip().split("\t")
        comment_dict["gyrA"][allele] = comment
with open(os.path.join(snakemake.params.mlst_dir, "db", "pubmlst", "ngstar", "NG_parC.tfa.comments")) as f:
    for line in f:
        allele, comment = line.rstrip().split("\t")
        comment_dict["parC"][allele] = comment
with open(os.path.join(snakemake.params.mlst_dir, "db", "pubmlst", "ngstar", "NG_23S.tfa.comments")) as f:
    for line in f:
        allele, comment = line.rstrip().split("\t")
        comment_dict["rna23S"][allele] = comment


outstring += "\t" + profile


allele = penA.split('(')[1].split(')')[0]
outstring += "\t" + allele
comment = []
for i in allele.split(','):
    comment.append(comment_dict["penA"][i.replace("~", "").replace("?", "")])
outstring += "\t" + ",".join(comment)

allele = mtrR.split('(')[1].split(')')[0]
outstring += "\t" + allele
comment = []
for i in allele.split(','):
    comment.append(comment_dict["mtrR"][i.replace("~", "").replace("?", "")])
outstring += "\t" + ",".join(comment)

allele = porB.split('(')[1].split(')')[0]
outstring += "\t" + allele
comment = []
for i in allele.split(','):
    comment.append(comment_dict["porB"][i.replace("~", "").replace("?", "")])
outstring += "\t" + ",".join(comment)

allele = ponA.split('(')[1].split(')')[0]
outstring += "\t" + allele
comment = []
for i in allele.split(','):
    comment.append(comment_dict["ponA"][i.replace("~", "").replace("?", "")])
outstring += "\t" + ",".join(comment)

allele = gyrA.split('(')[1].split(')')[0]
outstring += "\t" + allele
comment = []
for i in allele.split(','):
    comment.append(comment_dict["gyrA"][i.replace("~", "").replace("?", "")])
outstring += "\t" + ",".join(comment)

allele = parC.split('(')[1].split(')')[0]
outstring += "\t" + allele
comment = []
for i in allele.split(','):
    comment.append(comment_dict["parC"][i.replace("~", "").replace("?", "")])
outstring += "\t" + ",".join(comment)

allele = rna23S.split('(')[1].split(')')[0]
outstring += "\t" + allele
comment = []
for i in allele.split(','):
    comment.append(comment_dict["rna23S"][allele.replace("~", "").replace("?", "")])
outstring += "\t" + ",".join(comment)


with open(snakemake.input.ngmast) as f:
    contig, scheme, profile, porB, tbpB = f.readline().rstrip().split("\t")




outstring += "\t" + profile
outstring += "\t" + porB.split('(')[1].split(')')[0]
outstring += "\t" + tbpB.split('(')[1].split(')')[0]

rplf_dict = {}
with open(os.path.join(snakemake.params.mlst_dir, "db", "pubmlst", "rplf", "rplf.txt")) as f:
    for line in f:
        profile, rplf, species, comment = line.rstrip().split("\t")
        rplf_dict[profile] = [profile, species, comment]


with open(snakemake.input.rplf) as f:
    fasta, scheme, profile, rplf = f.readline().rstrip().split("\t")


if profile in rplf_dict:
    outstring += "\t" + "\t".join(rplf_dict[profile])
else:
    outstring += "\t" + profile + "\tmissing\tmissing"

with open(snakemake.input.ppng_cov) as f:
    f.readline()
    cov = f.readline().rstrip().split()[1]
    depth = f.readline().rstrip().split()[1]

outstring += "\t" + cov + "\t" + depth




with open(snakemake.output.tsv, 'w') as o:
    o.write("\t".join(headers) + "\n")
    o.write(outstring + "\n")
