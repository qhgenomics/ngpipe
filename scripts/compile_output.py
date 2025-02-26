import os


headers = ['Sample', 'NG-STAR', 'penA', 'mtrR', 'porB', 'ponA', 'gyrA', 'parC', '23S', 'NG-STAR_CC',
           'penA_mosaic_type', 'MLST', 'abcZ', 'adk', 'aroE', 'fumC', 'gdh', 'pdhC', 'pgm', 'NG-MAST', 'POR',
           'TBPB', 'rplF', 'rplF_species', 'rplf_depth', 'ppnG coverage', 'ppnG depth']

for i in snakemake.params.positions.split(','):
    headers.append('23S_bases_pos{}:a:t:g:c'.format(i))

outstring = snakemake.params.sample

with open(snakemake.input.pyngo) as f:
    header = f.readline().rstrip()
    body = f.readline().rstrip().split("\t")
    outstring += "\t" + "\t".join(body[1:])

rplf_dict = {}
with open(os.path.join(snakemake.params.mlst_dir, "rplf", "species.txt")) as f:
    for line in f:
        profile, species = line.rstrip().split("\t")
        rplf_dict[profile] = species


with open(snakemake.input.rplf) as f:
    f.readline()
    fasta, profile, rplf = f.readline().rstrip().split("\t")


if profile in rplf_dict:
    outstring += "\t" + profile + "\t" + rplf_dict[profile]
else:
    outstring += "\t" + profile + "\tmissing"

with open(snakemake.input.rplf_cov) as f:
    f.readline()
    cov = f.readline().rstrip().split()[1]
    depth = f.readline().rstrip().split()[1]

outstring += "\t" + depth


with open(snakemake.input.ppng_cov) as f:
    f.readline()
    cov = f.readline().rstrip().split()[1]
    depth = f.readline().rstrip().split()[1]

outstring += "\t" + cov + "\t" + depth

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
