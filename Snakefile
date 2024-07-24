configfile: workflow.source_path("config.yaml")
workdir: config["workdir"]


onsuccess:
    print("Workflow finished, no error")

onerror:
    print("An error occurred")





def get_samples(filename):
    with open(filename) as f:
        sample_list = []
        for line in f:
            if line.rstrip() != "":
                sample_list.append(line.rstrip())
    return sample_list

rule all:
    input:
        summary = "summary.tsv",
        mst = "mst.svg"



rule create_mst:
    input:
        summary = "summary.tsv"
    output:
        mst = "mst.svg"
    params:
        previous_runs = config["previous_runs"]
    script:
        "scripts/create_mst.py"


rule compile_output:
    input:
        summaries = expand("step6_output/{sample}_summary.tsv", sample=get_samples(config["samples"]))
    output:
        summary = "summary.tsv"
    run:
        first = True
        with open(output.summary, 'w') as o:
            for i in input.summaries:
                with open(i) as f:
                    if first:
                        first = False
                        o.write(f.readline())
                    else:
                        f.readline()
                    o.write(f.readline())








rule update_ngstar:
    params:
        mlst_dir = config["mlst_dir"],
        update_db = config["update_db"],
        db = "ngstar"
    output:
        log = "ngstar.log"
    script:
        "scripts/update_db.py"

rule update_ngmast:
    params:
        mlst_dir = config["mlst_dir"],
        update_db = config["update_db"],
        db = "ngmast"
    output:
        log = "ngmast.log"
    script:
        "scripts/update_db.py"


rule update_rplf:
    params:
        mlst_dir = config["mlst_dir"],
        update_db = config["update_db"],
        db = "rplf"
    output:
        log = "rplf.log"
    script:
        "scripts/update_db.py"



rule update_mlst:
    params:
        mlst_dir = config["mlst_dir"],
        update_db = config["update_db"],
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
        rplf_log = "rplf.log",
        mlst_log = "mlst.log"
    output:
        updated_db = "database.log"
    shell:
        "sed -i \"s|NG-MAST_||g\" {params.mlst_dir}/db/pubmlst/ngmast/* && "
        "sed -i \"s|'mtrR|mtrR|g\" {params.mlst_dir}/db/pubmlst/ngstar/* && "
        "sed -i \"s|'rplF|rplF|g\" {params.mlst_dir}/db/pubmlst/rplf/* && "
        "sed -i \"s|rplF_id|ST|g\" {params.mlst_dir}/db/pubmlst/rplf/rplf.txt &&"
        "sed -i \"s|genospecies|species|g\" {params.mlst_dir}/db/pubmlst/rplf/rplf.txt &&"
        "sed -i \"s|comments|CC|g\" {params.mlst_dir}/db/pubmlst/rplf/rplf.txt &&"
        "{params.mlst_dir}/scripts/mlst-make_blast_db && "
        "cat {input.ngstar_log} {input.ngmast_log} {input.mlst_log} {input.rplf_log} > {output.updated_db}"

rule qc:
    params:
        read_dir = config["read_dir"]
    output:
        qc1 = "step1_fastqc/{sample}_R1_fastqc.html",
        qc2 = "step1_fastqc/{sample}_R2_fastqc.html"
    threads: 24
    run:
        import subprocess
        if params.read_dir != "none":
            subprocess.Popen("fastqc {}/{}_R1.fastq.gz {}/{}_R2.fastq.gz -o step1_fastqc -t {}".format(
                params.read_dir, wildcards.sample, params.read_dir, wildcards.sample,threads), shell=True).wait()
        else:
            subprocess.Popen("mkdir -p step1_fastqc && touch {} && touch {}".format(output.qc1, output.qc2), shell=True).wait()

rule assemble_reads:
    params:
        read_dir = config["read_dir"],
        contig_dir = config["contig_dir"]
    output:
        scaffolds = "step2_assembly/metaspades_{sample}/scaffolds.fasta"
    threads: 24
    run:
        import subprocess
        if params.contig_dir != "none":
            subprocess.Popen(
                "mkdir -p step2_assembly/metaspades_{} && cp {}/{}.fasta {}".format(
                    wildcards.sample, params.contig_dir, wildcards.sample, output.scaffolds), shell=True).wait()
        else:
            try:
                subprocess.Popen(
                    "spades.py --meta -k 21,31,41,51,61,71,81,91,101,111 -o step2_assembly/metaspades_{} \
                    -1 {}/{}_R1.fastq.gz -2 {}/{}_R2.fastq.gz -t {}".format(
                    wildcards.sample, params.read_dir, wildcards.sample, params.read_dir, wildcards.sample, threads),
                shell=True).wait()
            except subprocess.CalledProcessError:
                pass
            if not os.path.exists(output.scaffolds):
                with open(output.scaffolds, 'w') as o:
                    o.write(">ntc\nnnnnnnnnnnnnnnnnnnnnnnnnnnnnnnn\n")



rule ng_typing:
    input:
        updated_db = "database.log",
        scaffolds = "step2_assembly/metaspades_{sample}/scaffolds.fasta"
    params:
        mlst_dir = config["mlst_dir"]
    output:
        mlst = "step3_typing/{sample}_mlst.tsv",
        ngmast = "step3_typing/{sample}_ngmast.tsv",
        ngstar = "step3_typing/{sample}_ngstar.tsv",
        rplf = "step3_typing/{sample}_rplf.tsv"
    threads: 24
    shell:
        "{params.mlst_dir}/bin/mlst --scheme mlst --threads 32 --quiet {input.scaffolds} > {output.mlst} & "
        "{params.mlst_dir}/bin/mlst --scheme ngmast --threads 32 --quiet {input.scaffolds} > {output.ngmast} & "
        "{params.mlst_dir}/bin/mlst --scheme rplf --threads 32 --quiet {input.scaffolds} > {output.rplf} & "
        "mv {params.mlst_dir}/db/pubmlst/ngmast/porB.tfa {params.mlst_dir}/db/pubmlst/ngmast/porB.tfa.hide && "
        "{params.mlst_dir}/scripts/mlst-make_blast_db && "
        "{params.mlst_dir}/bin/mlst --scheme ngstar --threads 32 --quiet {input.scaffolds} > {output.ngstar} && "
        "mv {params.mlst_dir}/db/pubmlst/ngmast/porB.tfa.hide {params.mlst_dir}/db/pubmlst/ngmast/porB.tfa && "
        "{params.mlst_dir}/scripts/mlst-make_blast_db"

rule abricate:
    input:
        scaffolds = "step2_assembly/metaspades_{sample}/scaffolds.fasta"
    output:
        abricate = "step4_abricate/{sample}_abricate.txt"
    threads: 24
    shell:
        "abricate {input.scaffolds} > {output.abricate}"


rule ppng_mapping:
    params:
        read_dir = config["read_dir"],
        contig_dir = config["contig_dir"],
        reference = config["ppng_fasta"]
    output:
        cov = "step5_cov/{sample}_ppng_coverage.txt",
        bam = "step5_cov/{sample}_ppng.bam"
    run:
        import subprocess
        if params.read_dir != "none":
            subprocess.Popen("minimap2 -ax sr {} {}/{}_R1.fastq.gz {}/{}_R2.fastq.gz | samtools view -bS - | "
                             "samtools sort -o {} && samtools depth -aa {} > {}".format(
                params.reference, params.read_dir, wildcards.sample, params.read_dir,
                wildcards.sample, output.bam, output.bam, output.cov), shell=True).wait()
        else:
            subprocess.Popen("minimap2 -ax asm5 {} {}/{}.fasta | samtools view -bS - | samtools sort -o {} && "
                             "samtools depth -aa {} > {}".format(params.reference, params.contig_dir,
                wildcards.sample, output.bam, output.bam, output.cov), shell=True).wait()


rule rplf_mapping:
    params:
        read_dir = config["read_dir"],
        contig_dir = config["contig_dir"],
        reference = config["rplf_fasta"]
    output:
        cov = "step5_cov/{sample}_rplf_coverage.txt",
        bam = "step5_cov/{sample}_rplf.bam"
    run:
        import subprocess
        if params.read_dir != "none":
            subprocess.Popen("minimap2 -ax sr {} {}/{}_R1.fastq.gz {}/{}_R2.fastq.gz | samtools view -bS - | "
                             "samtools sort -o {} && samtools depth -aa {} > {}".format(
                params.reference, params.read_dir, wildcards.sample, params.read_dir,
                wildcards.sample, output.bam, output.bam, output.cov), shell=True).wait()
        else:
            subprocess.Popen("minimap2 -ax asm5 {} {}/{}.fasta | samtools view -bS - | samtools sort -o {} && "
                             "samtools depth -aa {} > {}".format(params.reference, params.contig_dir,
                wildcards.sample, output.bam, output.bam, output.cov), shell=True).wait()

rule rrna_mapping:
    params:
        read_dir = config["read_dir"],
        contig_dir = config["contig_dir"],
        reference = config["rrna_fasta"]
    output:
        bam = "step5_cov/{sample}_rrna.bam"
    run:
        import subprocess
        if params.read_dir != "none":
            subprocess.Popen("minimap2 -ax sr {} {}/{}_R1.fastq.gz {}/{}_R2.fastq.gz | samtools view -bS - | "
                             "samtools sort -o {} && samtools index {}".format(
                params.reference, params.read_dir, wildcards.sample, params.read_dir,
                wildcards.sample, output.bam, output.bam), shell=True).wait()
        else:
            subprocess.Popen("minimap2 -ax asm5 {} {}/{}.fasta | samtools view -bS - | samtools sort -o {} && "
                             "samtools index {}".format(params.reference, params.contig_dir,
                              wildcards.sample, output.bam, output.bam), shell=True).wait()


rule rrna_alleles:
    params:
        reference = config["rrna_fasta"],
        position1 = config["position1"],
        position2 = config["position2"]
    input:
        bam = "step5_cov/{sample}_rrna.bam"
    output:
        tsv = "step5_cov/{sample}_rrna_alleles.tsv"
    run:
        import pysam
        samfile = pysam.AlignmentFile(input.bam, "rb")
        freqdict1 = {"a":0, "t":0, "c":0, "g":0}
        for pileupcolumn in samfile.pileup("23S", params.position1):
            if pileupcolumn.pos != params.position1:
                continue
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    try:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                    except IndexError:
                        pass
                freqdict1[base.lower()] += 1
        freqdict2 = {"a":0, "t":0, "c":0, "g":0}
        for pileupcolumn in samfile.pileup("23S", params.position2):
            if pileupcolumn.pos != params.position2:
                continue
            for pileupread in pileupcolumn.pileups:
                if not pileupread.is_del and not pileupread.is_refskip:
                    try:
                        base = pileupread.alignment.query_sequence[pileupread.query_position]
                    except IndexError:
                        pass
                freqdict2[base.lower()] += 1
        with open(output.tsv, 'w') as o:
            o.write("pos\ta\tt\tg\tc\n")
            o.write("{}\t{}\t{}\t{}\t{}\n".format(params.position1, freqdict1['a'], freqdict1['t'], freqdict1['g'], freqdict1['c']))
            o.write("{}\t{}\t{}\t{}\t{}\n".format(params.position2, freqdict2['a'], freqdict2['t'], freqdict2['g'], freqdict2['c']))


rule get_rplf_coverage:
    input:
        coverage = "step5_cov/{sample}_rplf_coverage.txt"
    output:
        stats = "step5_cov/{sample}_rplf.cov"
    run:
        with open(input.coverage) as f, open(output.stats, 'w') as o:
            total_depth, total_cov, total = 0, 0, 0
            for line in f:
                ref, pos, cov = line.split()
                total += 1
                total_depth += float(cov)
                if cov != '0':
                    total_cov += 1
            if total == 0:
                total = 1
            o.write("Bases: {}\nCoverage: {:.2%}\nDepth: {:.2f}".format(total, total_cov/total, total_depth/total))

rule get_ppng_coverage:
    input:
        coverage = "step5_cov/{sample}_ppng_coverage.txt"
    output:
        stats = "step5_cov/{sample}_ppng.cov"
    run:
        with open(input.coverage) as f, open(output.stats, 'w') as o:
            total_depth, total_cov, total = 0, 0, 0
            for line in f:
                ref, pos, cov = line.split()
                total += 1
                total_depth += float(cov)
                if cov != '0':
                    total_cov += 1
            if total == 0:
                total = 1
            o.write("Bases: {}\nCoverage: {:.2%}\nDepth: {:.2f}".format(total, total_cov/total, total_depth/total))



rule create_output:
    input:
        updated_db = "database.log",
        qc1 = "step1_fastqc/{sample}_R1_fastqc.html",
        qc2 = "step1_fastqc/{sample}_R2_fastqc.html",
        scaffolds = "step2_assembly/metaspades_{sample}/scaffolds.fasta",
        mlst = "step3_typing/{sample}_mlst.tsv",
        ngmast = "step3_typing/{sample}_ngmast.tsv",
        ngstar = "step3_typing/{sample}_ngstar.tsv",
        rplf = "step3_typing/{sample}_rplf.tsv",
        abricate = "step4_abricate/{sample}_abricate.txt",
        ppng_cov = "step5_cov/{sample}_ppng.cov",
        rplf_cov = "step5_cov/{sample}_rplf.cov",
        rrna_alleles = "step5_cov/{sample}_rrna_alleles.tsv"
    params:
        mlst_dir = config["mlst_dir"],
        sample = "{sample}",
        position1 = config["position1"],
        position2 = config["position2"]
    output:
        tsv = "step6_output/{sample}_summary.tsv"
    script:
        "scripts/compile_output.py"
