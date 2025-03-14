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
        results = "results.tsv",
        qc = "qclog.tsv",
        qc_summary = "summary.tsv",
        multi_qc_report = "multiqc_report/multiqc_report.html",
        mst = "mst.svg"



rule create_mst:
    input:
        summary = "results.tsv"
    output:
        mst = "mst.svg"
    params:
        previous_runs = config["previous_runs"],
        max_missing = config["max_missing_mst"]
    script:
        "scripts/create_mst.py"

rule multi_qc:
    input:
        qc = "qclog.tsv"
    output:
        multi_qc_report = "multiqc_report/multiqc_report.html"
    shell:
        "multiqc step1_qc/ -e samtools -o multiqc_report -f"


rule compile_output:
    input:
        summaries = expand("step6_output/{sample}_summary.tsv", sample=get_samples(config["samples"]))
    output:
        results = "results.tsv",
        qc = "qclog.tsv",
        summary = "summary.tsv"
    run:
        first = True
        with open(output.results, 'w') as o, open(output.qc, 'w') as qc:
            for i in input.summaries:
                with open(i) as f:
                    if first:
                        first = False
                        header = f.readline().rstrip().split("\t")
                        for num, j in enumerate(header):
                            if j.startswith("23S_bases_pos"):
                                qcsplit = num
                        qcsplit += 1
                        freq_dict = {}
                        for i in header[1:qcsplit]:
                            freq_dict[i] = {}
                        o.write("\t".join(header[:qcsplit]) + "\n")
                        qc.write("Sample\t" + "\t".join(header[qcsplit:]) + "\n")
                    else:
                        f.readline()
                    body = f.readline().rstrip().split("\t")
                    o.write("\t".join(body[:qcsplit]) + "\n")
                    for cat, j in zip(header[1:qcsplit], body[1:qcsplit]):
                        if not j in freq_dict[cat]:
                            freq_dict[cat][j] = 0
                        freq_dict[cat][j] += 1
                    qc.write("{}\t".format(body[0]) + "\t".join(body[qcsplit:]) + "\n")
        with open(output.summary,'w') as o:
            outstring = ""
            cats = {}
            maxcats = 0
            for i in header[1:qcsplit]:
                outstring += i + "_category\t" + i + "_count\t" + i + "_percent\t"
                thelist = list(freq_dict[i])
                thelist.sort(key=lambda x: freq_dict[i][x], reverse=True)
                cats[i] = thelist
                if len(thelist) > maxcats:
                    maxcats = len(thelist)
            outstring = outstring[:-1] + "\n"
            for i in range(maxcats):
                for j in header[1:qcsplit]:
                    if i < len(cats[j]):
                        outstring += cats[j][i] + "\t" + str(freq_dict[j][cats[j][i]]) + "\t{:.2%}\t".format(freq_dict[j][cats[j][i]]/len(get_samples(config["samples"])))
                    else:
                        outstring += "\t\t\t"
                outstring = outstring[:-1] + "\n"
            o.write(outstring)






rule update_rplf:
    params:
        mlst_dir = config["mlst_db"],
        update_db = config["update_db"],
        db = "rplf"
    output:
        log = "rplf.log"
    script:
        "scripts/update_db.py"



rule update_mlst:
    params:
        mlst_dir = config["mlst_db"],
        update_db = config["update_db"],
        db = "mlst"
    output:
        log = "mlst.log",
    script:
        "scripts/update_db.py"


rule update_ngstar:
    params:
        mlst_dir = config["mlst_db"],
        update_db = config["update_db"],
        db = "ngstar"
    output:
        log = "ngstar.log",
    script:
        "scripts/update_db.py"

rule update_db:
    input:
        mlst = "mlst.log",
        rplf = "rplf.log",
        ngstar = "ngstar.log"
    params:
        mlst_db = config["mlst_db"],
        update_db = config["update_db"],
        ngstar_cc = workflow.source_path("data/240513_NGSTAR_CC_updated.csv")
    output:
        updated_db = "database.log",
    run:
        from datetime import datetime
        download_date_file = os.path.join(params.mlst_db, "download_date.log")
        today = datetime.today().strftime('%Y-%m-%d')
        download_date = None
        if os.path.exists(download_date_file):
            with open(download_date_file) as f:
                download_date = f.readline().split(":")[0]
        if params.update_db and download_date == today and os.path.exists(params.mlst_db):
            log_text = "{}: Database already downloaded today, skipping download.".format(today)
        elif params.update_db:
            if download_date is None:
                log_text = "{}: no database found. Downloading database".format(today)
            else:
                log_text = "{}: found database updated {}, updating.".format(today, download_date)
            shell("pyngoST.py -d -n {} -cc {}".format(params.mlst_db, params.ngstar_cc))
            with open(download_date_file, 'w') as f:
                f.write("{}: database downloaded.".format(today))
        elif download_date is None or not os.path.exists(params.mlst_db):
            raise ValueError("Database file missing, please point to pyngost database or set update_db to True.")
        else:
            log_text = "Database update not requested: using database downloaded on {}.".format(download_date)
        shell("claMLST create --force {params.mlst_db}/pymlst_rplf {params.mlst_db}/rplf/rplf.txt {params.mlst_db}/rplf/rplF.tfa")
        shell("claMLST create --force {params.mlst_db}/pymlst_mlst {params.mlst_db}/mlst/mlst.txt {params.mlst_db}/mlst/abcZ.tfa "
              "{params.mlst_db}/mlst/adk.tfa {params.mlst_db}/mlst/aroE.tfa {params.mlst_db}/mlst/fumC.tfa "
              "{params.mlst_db}/mlst/gdh.tfa {params.mlst_db}/mlst/pdhC.tfa {params.mlst_db}/mlst/pgm.tfa ")
        with open(output.updated_db, 'w') as o, open(input.mlst) as mlst_log, open(input.rplf) as rplf_log:
            o.write("mlst updated: " + mlst_log.readline())
            o.write("rplf updated: " + rplf_log.readline())
            o.write(log_text)





rule qc:
    params:
        read_dir = config["read_dir"],
        quast_fasta = workflow.source_path("data/ref.fasta"),
        quast_gff = workflow.source_path("data/ref.gff3")
    input:
        scaffolds = "step2_assembly/metaspades_{sample}/scaffolds.fasta"
    output:
        qc1 = "step1_qc/{sample}_R1_fastqc.html",
        qc2 = "step1_qc/{sample}_R2_fastqc.html",
        quast = "step1_qc/{sample}_quast/report.txt"
    threads: 24
    run:
        import subprocess
        read1 = "{}/{}_R1.fastq.gz".format(params.read_dir, wildcards.sample)
        read2 = "{}/{}_R2.fastq.gz".format(params.read_dir, wildcards.sample)
        if params.read_dir != "none":
            subprocess.Popen("fastqc {} {} -o step1_qc -t {}".format(read1, read2, threads), shell=True).wait()
            subprocess.Popen("quast.py {} -r {} -g {} -L -o step1_qc/{}_quast/ -1 {} -2 {} --threads {}".format(
                input.scaffolds, params.quast_fasta, params.quast_gff, wildcards.sample, read1, read2, threads), shell=True).wait()
        else:
            subprocess.Popen("mkdir -p step1_fastqc && touch {} && touch {}".format(output.qc1, output.qc2), shell=True).wait()
            subprocess.Popen("quast.py {} -r {} -g {} -L -o step1_qc/{}_quast/ --threads {}".format(
                input.scaffolds, params.quast_fasta, params.quast_gff, wildcards.sample, threads), shell=True).wait()

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
        scaffolds = "step2_assembly/metaspades_{sample}/scaffolds.fasta",
    params:
        mlst_db = config["mlst_db"]
    output:
        pyngo_out = "step3_typing/{sample}_pyngo.tsv",
        mlst= "step3_typing/{sample}_mlst.tsv",
        rplf= "step3_typing/{sample}_rplf.tsv"
    threads: 24
    shell:
        "pyngoST.py -i {input.scaffolds} -s NG-STAR,MLST,NG-MAST -c -m -b -a -o {wildcards.sample}_pyngo.tsv -q step3_typing -t {threads} -p {params.mlst_db} && "
        "claMLST search -o {output.mlst} {params.mlst_db}/pymlst_mlst {input.scaffolds} && "
        "claMLST search -o {output.rplf} {params.mlst_db}/pymlst_rplf {input.scaffolds}"

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
        reference = workflow.source_path("data/ppng_target.fasta")
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
        reference = workflow.source_path("data/rplf_target.fasta")
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
        reference = workflow.source_path("data/rrna_target.fasta")
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
        reference = workflow.source_path("data/rrna_target.fasta"),
        positions = config["positions"]
    input:
        bam = "step5_cov/{sample}_rrna.bam"
    output:
        tsv = "step5_cov/{sample}_rrna_alleles.tsv"
    run:
        import pysam
        samfile = pysam.AlignmentFile(input.bam, "rb")
        freqdicts = []
        for position in params.positions.split(','):
            freqdict = {"a":0, "t":0, "c":0, "g":0}
            for pileupcolumn in samfile.pileup("23S", int(position)):
                if pileupcolumn.pos != int(position):
                    continue
                for pileupread in pileupcolumn.pileups:
                    if not pileupread.is_del and not pileupread.is_refskip:
                        try:
                            base = pileupread.alignment.query_sequence[pileupread.query_position]
                        except IndexError:
                            pass
                    freqdict[base.lower()] += 1
            freqdicts.append(freqdict)
        with open(output.tsv, 'w') as o:
            o.write("pos\ta\tt\tg\tc\n")
            for freqdict, position in zip(freqdicts, params.positions.split(',')):
                o.write("{}\t{}\t{}\t{}\t{}\n".format(position, freqdict['a'], freqdict['t'], freqdict['g'], freqdict['c']))
          

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


rule get_target_coverage:
    input:
        pyngo = "step3_typing/{sample}_pyngo.tsv",
    params:
        read_dir = config["read_dir"],
        contig_dir = config["contig_dir"],
        mlst_db = config["mlst_db"]
    output:
        mlst_cov = "step5_cov/{sample}.mlst.cov",
        ngstar_cov = "step5_cov/{sample}.ngstar.cov",
        ngmast_cov = "step5_cov/{sample}.ngmast.cov",

    script:
        "scripts/get_loci_coverage.py"

rule create_output:
    input:
        updated_db = "database.log",
        qc1 = "step1_qc/{sample}_R1_fastqc.html",
        qc2 = "step1_qc/{sample}_R2_fastqc.html",
        scaffolds = "step2_assembly/metaspades_{sample}/scaffolds.fasta",
        pyngo = "step3_typing/{sample}_pyngo.tsv",
        rplf = "step3_typing/{sample}_rplf.tsv",
        abricate = "step4_abricate/{sample}_abricate.txt",
        ppng_cov = "step5_cov/{sample}_ppng.cov",
        rplf_cov = "step5_cov/{sample}_rplf.cov",
        rrna_alleles = "step5_cov/{sample}_rrna_alleles.tsv",
        mlst_cov = "step5_cov/{sample}.mlst.cov",
        ngstar_cov = "step5_cov/{sample}.ngstar.cov",
        ngmast_cov = "step5_cov/{sample}.ngmast.cov",
    params:
        mlst_dir = config["mlst_db"],
        sample = "{sample}",
        positions = config["positions"]
    output:
        tsv = "step6_output/{sample}_summary.tsv"
    script:
        "scripts/compile_output.py"
