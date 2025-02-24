import subprocess
import os



def create_coverage_files(input_file, scheme, mlst_dir, outfile):
    with (open(input_file) as f, open("step5_cov/{}.{}.fasta".format(snakemake.wildcards.sample, scheme), "w") as o):
        header = f.readline().rstrip().split("\t")
        splitline = f.readline().rstrip().split("\t")
        contig, profile = splitline[:2]
        alleles = splitline[2:]
        for allele, gene in zip(alleles, header[2:]):
            if ',' in allele:
                allele = allele.split(',')[0]
            allele = allele.replace("?", "").replace("~", "")
            getseq = False
            if allele == "" or allele == "-" or allele == "new":
                allele = "1"
            with open("{}/{}/{}.tfa".format(mlst_dir, scheme, gene)) as f:
                for line in f:
                    if line.startswith(">"):
                        if getseq:
                            getseq = False
                            break
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

create_coverage_files(snakemake.input.mlst, "mlst", snakemake.params.mlst_dir, snakemake.output.mlst_cov )
create_coverage_files(snakemake.input.ngstar, "ngstar", snakemake.params.mlst_dir, snakemake.output.ngstar_cov)
create_coverage_files(snakemake.input.ngmast, "ngmast", snakemake.params.mlst_dir, snakemake.output.ngmast_cov)
