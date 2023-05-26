configfile: workflow.source_path("config.yaml")
workdir: config["workdir"]


onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")


rule all:
    input:
        expand("step7_output/{sample}_summary.tsv", sample=config["sample"])


rule update_ngstar:
    params:
        mlst_dir = config["mlst_dir"],
        db = "ngstar"
    output:
        log = "ngstar.log"
    script:
        "scripts/update_db.py"

rule update_ngmast:
    params:
        mlst_dir = config["mlst_dir"],
        db = "ngmast"
    output:
        log = "ngmast.log"
    script:
        "scripts/update_db.py"

rule update_mlst:
    params:
        mlst_dir = config["mlst_dir"],
        db = "mlst"
    output:
        log = "mlst.log",
    script:
        "scripts/update_db.py"


rule update_db:
    params:
        mlst_dir = config["mlst_dir"]
    input:
        ngstar_log = "ngstar.log",
        ngmast_log = "ngmast.log",
        mlst_log = "mlst.log"
    output:
        updated_db = "database.log"
    shell:
        "sed -i \"s|NG-MAST_||g\" {params.mlst_dir}/db/pubmlst/ngmast/* && "
        "sed -i \"s|'mtrR|mtrR|g\" {params.mlst_dir}/db/pubmlst/ngstar/* && "
        "{params.mlst_dir}/scripts/mlst-make_blast_db && "
        "cat {input.ngstar_log} {input.ngmast_log} {input.mlst_log} > {output.updated_db}"

rule qc:
    params:
        read_dir = config["read_dir"],
        sample = config["sample"]
    output:
        qc1 = "step1_fastqc/{sample}_R1_fastqc.html",
        qc2 = "step1_fastqc/{sample}_R2_fastqc.html"
    threads: 24
    shell:
        "fastqc {params.read_dir}/{params.sample}_R1.fastq.gz {params.read_dir}/{params.sample}_R2.fastq.gz -o step1_fastqc -t {threads}"

rule assemble_reads:
    params:
        read_dir = config["read_dir"],
        sample = config["sample"],
    output:
        scaffolds = "step2_assembly_metaspades_{sample}/scaffolds.fasta"
    threads: 24
    shell:
        "spades.py --meta -k 21,31,41,51,61,71,81,91,101,111 -o step2_assembly_metaspades_{params.sample} -1 {params.read_dir}/{params.sample}_R1.fastq.gz -2 {params.read_dir}/{params.sample}_R2.fastq.gz -t {threads}"

rule ng_typing:
    input:
        updated_db = "database.log",
        scaffolds = "step2_assembly_metaspades_{sample}/scaffolds.fasta"
    params:
        mlst_dir = config["mlst_dir"]
    output:
        mlst = "step3_typing/{sample}_mlst.tsv",
        ngmast = "step3_typing/{sample}_ngmast.tsv",
        ngstar = "step3_typing/{sample}_ngstar.tsv"
    threads: 24
    shell:
        "{params.mlst_dir}/bin/mlst --scheme mlst --threads 32 --quiet {input.scaffolds} > {output.mlst} & "
        "{params.mlst_dir}/bin/mlst --scheme ngmast --threads 32 --quiet {input.scaffolds} > {output.ngmast} & "
        "{params.mlst_dir}/bin/mlst --scheme ngstar --threads 32 --quiet {input.scaffolds} > {output.ngstar}"

rule abricate:
    input:
        scaffolds = "step2_assembly_metaspades_{sample}/scaffolds.fasta"
    output:
        abricate = "step4_abricate/{sample}_abricate.txt"
    threads: 24
    shell:
        "abricate {input.scaffolds} > {output.abricate}"


rule ppng_mapping:
    params:
        read_dir = config["read_dir"],
        sample = config["sample"],
        reference = workflow.source_path("data/ppng_target.fasta")
    output:
        cov = "step5_ppng/{sample}_ppng_coverage.txt",
        bam = "step5_ppng/{sample}_ppng.bam"
    shell:
        "minimap2 -ax sr {params.reference} {params.read_dir}/{params.sample}_R1.fastq.gz {params.read_dir}/{params.sample}_R2.fastq.gz | samtools view -bS - | samtools sort -o {output.bam} && "
        "samtools depth -aa {output.bam} > {output.cov}"

rule rplf_mapping:
    params:
        read_dir = config["read_dir"],
        sample = config["sample"],
        reference = workflow.source_path("data/rplf_target.fasta")
    output:
        cov = "step6_rplf/{sample}_rplf_coverage.txt",
        bam = "step6_rplf/{sample}_rplf.bam"
    shell:
        "minimap2 -ax sr {params.reference} {params.read_dir}/{params.sample}_R1.fastq.gz {params.read_dir}/{params.sample}_R2.fastq.gz | samtools view -bS - | samtools sort -o {output.bam} && "
        "samtools depth -aa {output.bam} > {output.cov}"


rule get_gene_coverage:
    input: "step5_ppng/{sample}_ppng_coverage.txt", "step6_rplf/{sample}_rplf_coverage.txt"
    output: "step5_ppng/{sample}_ppng.cov","step6_rplf/{sample}_rplf.cov"
    run:
        for infile, outfile in zip(input, output):
            with open(infile) as f, open(outfile, 'w') as o:
                total_depth, total_cov, total = 0, 0, 0
                for line in f:
                    ref, pos, cov = line.split()
                    total += 1
                    total_depth += float(cov)
                    if cov != '0':
                        total_cov += 1
                if total == 0:
                    total = 1
                o.write("Bases: {}\nCoverage: {:.2%}\nDepth: {:.2}".format(total, total_cov/total, total_depth/total))



rule create_output:
    input:
        updated_db = "database.log",
        qc1 = "step1_fastqc/{sample}_R1_fastqc.html",
        qc2 = "step1_fastqc/{sample}_R2_fastqc.html",
        scaffolds = "step2_assembly_metaspades_{sample}/scaffolds.fasta",
        mlst = "step3_typing/{sample}_mlst.tsv",
        ngmast = "step3_typing/{sample}_ngmast.tsv",
        ngstar = "step3_typing/{sample}_ngstar.tsv",
        abricate = "step4_abricate/{sample}_abricate.txt",
        ppng_cov = "step5_ppng/{sample}_ppng.cov",
        rplf_cov = "step6_rplf/{sample}_rplf.cov"
    params:
        mlst_dir = config["mlst_dir"]
    output:
        tsv = "step7_output/{sample}_summary.tsv"
    script:
        "scripts/compile_output.py"