import subprocess
import os



def create_coverage_files(input_file, scheme, mlst_dir, outfile):
    with (open(input_file) as f, open("step5_cov/{}.{}.fasta".format(snakemake.wildcards.sample, scheme), "w") as o):
        splitline = f.readline().rstrip().split("\t")
        contig, scheme, profile = splitline[:3]
        alleles = splitline[3:]
        for i in alleles:
            gene = i.split('(')[0]
            allele = i.split('(')[1].split(')')[0]
            if ',' in allele:
                allele = allele.split(',')[0]
            allele.replace("?", "").replace("~", "")
            getseq = False
            with open("{}/db/pubmlst/{}/{}.tfa".format(mlst_dir.replace("bin/mlst", ""), scheme, gene)) as f:
                for line in f:
                    if line.startswith(">"):
                        if getseq:
                            getseq = False
                            break
                        elif allele == '-':
                            o.write(line)
                            getseq = True
                        elif line.rstrip() == ">{}_{}".format(gene, allele):
                            o.write(line)
                            getseq = True
                    elif getseq:
                        o.write(line)
    if snakemake.params.read_dir != "none":
        subprocess.Popen("minimap2 -ax sr step5_cov/{sample}.{scheme}.fasta {read_dir}/{sample}_R1.fastq.gz {read_dir}/{sample}_R2.fastq.gz"
                         " | samtools view -bS - | samtools sort -o step5_cov/{sample}.{scheme}.bam && "
                         " samtools depth -aa step5_cov/{sample}.{scheme}.bam > {cov}".format(
            sample=snakemake.wildcards.sample, scheme=scheme, read_dir=snakemake.params.read_dir, cov=outfile), shell=True).wait()
    else:
        subprocess.Popen("minimap2 -ax asm5 step5_cov/{sample}.{scheme}.fasta {contig_dir}/{sample}.fasta"
                         " | samtools view -bS - | samtools sort -o step5_cov/{sample}.{scheme}.bam && "
                         " samtools depth -aa step5_cov/{sample}.{scheme}.bam > {cov}".format(
            sample=snakemake.wildcards.sample, scheme=scheme, contig_dir=snakemake.params.contig_dir, cov=outfile), shell=True).wait()

create_coverage_files(snakemake.input.mlst, "mlst", snakemake.input.mlst_dir, snakemake.output.mlst_cov )
create_coverage_files(snakemake.input.mlst, "ngstar", snakemake.input.starmlst_dir, snakemake.output.ngstar_cov)
create_coverage_files(snakemake.input.mlst, "ngmast", snakemake.input.mlst_dir, snakemake.output.ngmast_cov)

