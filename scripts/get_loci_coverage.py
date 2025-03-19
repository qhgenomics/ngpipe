import subprocess
import os



def create_coverage_files(input_file, scheme, mlst_dir, outfile):
    with (open(input_file) as f, open("step5_cov/{}.{}.fasta".format(snakemake.wildcards.sample, scheme), "w") as o):
        header = f.readline().rstrip().split("\t")
        splitline = f.readline().rstrip().split("\t")
        alleles = splitline[2:]
        for allele, gene in zip(alleles, header[2:]):
            if '-' in allele:
                allele = allele.split('-')[0]
            if '|' in allele:
                allele = allele.split('|')[0]
            if "_" in allele:
                allele = allele.split("_")[0]
            if allele == "" or allele == "-" or allele == "new":
                allele = "1"
                if gene == "penA":
                    allele += '.001'
            if gene in ['penA', 'mtrR', 'porB', 'ponA', 'gyrA', 'parC', '23S']:
                gene_scheme = "ngstar"
            elif gene in ['abcZ', 'adk', 'aroE', 'fumC', 'gdh', 'pdhC', 'pgm']:
                gene_scheme = "mlst"
            elif gene in ['POR', 'TBPB']:
                gene_scheme = "ngmast"
            else:
                gene_scheme = None
            if gene_scheme != scheme:
                continue
            if scheme == "ngstar" and not '.' in allele:
                allele += '.0'
            if gene == "penA" and allele.endswith('0'):
                allele = allele[:-1]
            with open("{}/{}.fas".format(mlst_dir, gene)) as f:
                getseq = False
                for line in f:
                    if line.startswith(">"):
                        if getseq:
                            break
                        elif line.rstrip() == ">{}_{}".format(gene, allele):
                            o.write(line)
                            getseq = True
                    elif getseq:
                        o.write(line)
            if not getseq:
                raise Exception("{}_{} not found in {}/{}.fas".format(gene, allele, mlst_dir, gene))
    if snakemake.params.read_dir != "none":
        subprocess.Popen("minimap2 -a -k22 -w11 --sr --frag=yes -A2 -B8 -O12,32 -E2,1 -b0 -r100 -p.5 -N20 -f1000,5000 "
                     " -n1 -m20 -s40 -g10 -2K50m --heap-sort=yes --secondary=no step5_cov/{sample}.{scheme}.fasta "
                     "{read_dir}/{sample}_R1.fastq.gz {read_dir}/{sample}_R2.fastq.gz"
                     " | samtools view -bS - | samtools sort -o step5_cov/{sample}.{scheme}.bam && "
                     " samtools depth -aa step5_cov/{sample}.{scheme}.bam > {cov}".format(
        sample=snakemake.wildcards.sample, scheme=scheme, read_dir=snakemake.params.read_dir, cov=outfile), shell=True).wait()
    else:
         subprocess.Popen("minimap2 -a -k22 -w11 --sr --frag=yes -A2 -B8 -O12,32 -E2,1 -b0 -r100 -p.5 -N20 -f1000,5000 "
                          " -n1 -m20 -s40 -g10 -2K50m --heap-sort=yes --secondary=no step5_cov/{sample}.{scheme}.fasta "
                          "{contig_dir}/{sample}.fasta | samtools view -bS - | samtools sort -o step5_cov/{sample}.{scheme}.bam && "
                          " samtools depth -aa step5_cov/{sample}.{scheme}.bam > {cov}".format(
             sample=snakemake.wildcards.sample, scheme=scheme, contig_dir=snakemake.params.contig_dir, cov=outfile), shell=True).wait()

create_coverage_files(snakemake.input.pyngo, "mlst", snakemake.params.mlst_db, snakemake.output.mlst_cov )
create_coverage_files(snakemake.input.pyngo, "ngstar", snakemake.params.mlst_db, snakemake.output.ngstar_cov)
create_coverage_files(snakemake.input.pyngo, "ngmast", snakemake.params.mlst_db, snakemake.output.ngmast_cov)
