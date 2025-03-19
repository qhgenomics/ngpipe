import os


headers = ['Sample', 'NG-STAR', 'penA', 'mtrR', 'porB', 'ponA', 'gyrA', 'parC', '23S', 'NG-STAR_CC',
           'penA_mosaic_type', 'MLST', 'abcZ', 'adk', 'aroE', 'fumC', 'gdh', 'pdhC', 'pgm', 'NG-MAST', 'POR',
           'TBPB', 'rplF', 'rplF_species', 'rplf_coverage', 'rplf_depth', 'ppnG_present', 'ppnG coverage', 'ppnG depth',
           'ngstar_comment', 'penA_comment', 'mtrR_comment', 'porB_comment', 'ponA_comment', 'gyrA_comment', 'parC_comment', '23S_comment']

for i in snakemake.params.positions.split(','):
    headers.append('23S_bases_pos_{}:a:t:g:c'.format(i))

outstring = snakemake.params.sample

with open(snakemake.input.pyngo) as f:
    header = f.readline().rstrip()
    body = f.readline().rstrip().split("\t")
    ngstar_alleles = body[2:9]
    ngstar_st = body[1]
    outstring += "\t" + "\t".join(body[1:])

rplf_dict = {}
with open(os.path.join(snakemake.params.mlst_dir, "rplf", "species.txt")) as f:
    for line in f:
        profile, species = line.rstrip().split("\t")
        rplf_dict[profile] = species


with open(snakemake.input.rplf) as f:
    f.readline()
    fasta, profile, rplf = f.readline().rstrip("\n").split("\t")
    if profile == "":
        profile = '-'


if profile in rplf_dict:
    outstring += "\t" + profile + "\t" + rplf_dict[profile]
else:
    outstring += "\t" + profile + "\tmissing"

with open(snakemake.input.rplf_cov) as f:
    f.readline()
    cov = f.readline().rstrip().split()[1]
    depth = f.readline().rstrip().split()[1]

outstring += "\t" + cov + "\t" + depth


with open(snakemake.input.ppng_cov) as f:
    f.readline()
    cov = f.readline().rstrip().split()[1]
    depth = f.readline().rstrip().split()[1]

if (snakemake.params.read_dir == "none" or float(depth) >= 5) and float(cov[:-1]) >= 90:
    outstring += "\tTrue"
else:
    outstring += "\tFalse"

outstring += "\t" + cov + "\t" + depth

comment_file = os.path.join(snakemake.params.mlst_dir, "ngstar", "ngstar_profile.comments")
with open(comment_file) as f:
    comment = " not found"
    for line in f:
        st_file, iso_comment = line.rstrip().split("\t")[0:2]
        if st_file == ngstar_st:
            comment = iso_comment
            break
    outstring += "\t{}".format(comment)

for i, allele in zip(['penA', 'mtrR', 'porB', 'ponA', 'gyrA', 'parC', '23S'], ngstar_alleles):
    comment_file = os.path.join(snakemake.params.mlst_dir, "ngstar", i + ".comments")
    with open(comment_file) as f:
        f.readline()
        comment = "not found"
        for line in f:
            allele_file = line.split("\t")[0]
            iso_comment = line.split("\t")[1]
            if allele == allele_file:
                comment = iso_comment
        outstring += "\t{}".format(comment)


with open(snakemake.input.rrna_alleles) as f:
    f.readline()
    for line in f:
        pos, a, t, g, c = line.split()
        outstring += "\t{}:{}:{}:{}".format(a,t,g,c)

with open(snakemake.input.mlst_cov) as f:
    total_depth, total_cov, total = 0, 0, 0
    lastref = None
    for line in f:
        ref, pos, cov = line.split()
        if not lastref is None and ref != lastref:
            headers.append("{}_cov".format(lastref))
            headers.append("{}_depth".format(lastref))
            outstring += "\t{:.1%}\t{:.2f}".format(total_cov/total, total_depth/total)
            total_depth, total_cov, total = 0, 0, 0
        lastref = ref
        total += 1
        total_depth += float(cov)
        if cov != '0':
            total_cov += 1
    headers.append("{}_cov".format(lastref))
    headers.append("{}_depth".format(lastref))
    outstring += "\t{:.1%}\t{:.2f}".format(total_cov / total, total_depth / total)
    total_depth, total_cov, total = 0, 0, 0


with open(snakemake.input.ngstar_cov) as f:
    total_depth, total_cov, total = 0, 0, 0
    lastref = None
    for line in f:
        ref, pos, cov = line.split()
        if not lastref is None and ref != lastref:
            headers.append("{}_cov".format(lastref))
            headers.append("{}_depth".format(lastref))
            outstring += "\t{:.1%}\t{:.2f}".format(total_cov/total, total_depth/total)
            total_depth, total_cov, total = 0, 0, 0
        lastref = ref
        total += 1
        total_depth += float(cov)
        if cov != '0':
            total_cov += 1
    headers.append("{}_cov".format(lastref))
    headers.append("{}_depth".format(lastref))
    outstring += "\t{:.1%}\t{:.2f}".format(total_cov / total, total_depth / total)
    total_depth, total_cov, total = 0, 0, 0

with open(snakemake.input.ngmast_cov) as f:
    total_depth, total_cov, total = 0, 0, 0
    lastref = None
    for line in f:
        ref, pos, cov = line.split()
        if not lastref is None and ref != lastref:
            headers.append("{}_cov".format(lastref))
            headers.append("{}_depth".format(lastref))
            outstring += "\t{:.1%}\t{:.2f}".format(total_cov/total, total_depth/total)
            total_depth, total_cov, total = 0, 0, 0
        lastref = ref
        total += 1
        total_depth += float(cov)
        if cov != '0':
            total_cov += 1
    headers.append("{}_cov".format(lastref))
    headers.append("{}_depth".format(lastref))
    outstring += "\t{:.1%}\t{:.2f}".format(total_cov / total, total_depth / total)
    total_depth, total_cov, total = 0, 0, 0



with open(snakemake.output.tsv, 'w') as o:
    o.write("\t".join(headers) + "\n")
    o.write(outstring + "\n")
