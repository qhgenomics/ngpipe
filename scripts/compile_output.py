import os


headers = ['Sample', 'NG-STAR', 'penA', 'mtrR', 'porB', 'ponA', 'gyrA', 'parC', '23S', 'NG-STAR_CC',
           'penA_mosaic_type', 'MLST', 'abcZ', 'adk', 'aroE', 'fumC', 'gdh', 'pdhC', 'pgm', 'NG-MAST', 'POR',
           'TBPB', 'rplF', 'rplF_species', 'rplf_coverage', 'rplf_depth', 'ppnG_present', 'ppnG coverage', 'ppnG depth',
           'ngstar_comment', 'penA_comment', 'mtrR_comment', 'porB_comment', 'ponA_comment', 'gyrA_comment', 'parC_comment', '23S_comment']

for i in snakemake.params.positions.split(','):
    headers.append('23S_bases_pos_{}:a:t:g:c'.format(i))

outstring = snakemake.params.sample


pymlst_dict = {}
with open("step3_typing/{}_pymlst_mlst.tsv".format(snakemake.params.sample)) as f:
    header = f.readline().rstrip().split("\t")[2:]
    body = f.readline().rstrip().split("\t")[2:]
    for i, j in zip(header, body):
        if not ';' in j and not j == "new" and not j == '':
            pymlst_dict[i] = j

with open("step3_typing/{}_pymlst_ngstar.tsv".format(snakemake.params.sample)) as f:
    header = f.readline().rstrip().split("\t")[2:]
    body = f.readline().rstrip().split("\t")[2:]
    for i, j in zip(header, body):
        if not ';' in j and not j == "new" and not j == '':
            pymlst_dict[i.split('-')[0]] = j

with open("step3_typing/{}_pymlst_ngmast.tsv".format(snakemake.params.sample)) as f:
    header = f.readline().rstrip().split("\t")[2:]
    body = f.readline().rstrip().split("\t")[2:]
    for i, j in zip(header, body):
        if not ';' in j and not j == "new" and not j == '':
            pymlst_dict[i] = j

with open("step3_typing/{}_pymlst_rplf.tsv".format(snakemake.params.sample)) as f:
    header = f.readline().rstrip().split("\t")[2:]
    body = f.readline().rstrip().split("\t")[2:]
    for i, j in zip(header, body):
        if not ';' in j and not j == "new" and not j == '':
            pymlst_dict[i] = j



with open(snakemake.input.pyngo) as f:
    header = f.readline().rstrip().split("\t")
    body = f.readline().rstrip().split("\t")
    new_body = []
    for i, j in zip(header, body):
        if j.endswith("-1") and i in pymlst_dict:
            new_body.append(pymlst_dict[i])
        else:
            new_body.append(j)
    if body[2:9] != new_body[2:9]:
        new_st = '-'
        with open(os.path.join(snakemake.params.mlst_dir,  "NGSTAR_profiles.tab")) as profile:
            for line in profile:
                if line.split()[1:8] == new_body[2:9]:
                    new_st = line.split()[0]
                    break
        new_body[1] = new_st
        new_cc = '-'
        with open(snakemake.params.ngstar_cc) as cc:
            for line in cc:
                if line.split(',')[0] == new_st:
                    new_cc = line.rstrip().split(',')[1]
                    break
        new_body[9] = new_cc
    if body[12:19] != new_body[12:19]:
        new_st = '-'
        with open(os.path.join(snakemake.params.mlst_dir, "MLST_profiles.tab")) as profile:
            for line in profile:
                if line.split()[1:8] == new_body[12:19]:
                    new_st = line.split()[0]
                    break
        new_body[11] = new_st
    if body[20:22] != new_body[20:22]:
        new_st = '-'
        with open(os.path.join(snakemake.params.mlst_dir, "NGMAST_profiles.tab")) as profile:
            for line in profile:
                if line.split()[1:3] == new_body[20:22]:
                    new_st = line.split()[0]
                    break
        new_body[19] = new_st
    body = new_body
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
