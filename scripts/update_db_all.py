from datetime import datetime
import os
import subprocess

if not os.path.exists(os.path.join(snakemake.params.mlst_db, "pymlst")):
    os.makedirs(os.path.join(snakemake.params.mlst_db, "pymlst"))

download_date_file = os.path.join(snakemake.params.mlst_db, "download_date.log")
today = datetime.today().strftime('%Y-%m-%d')
download_date = None
if os.path.exists(download_date_file):
    with open(download_date_file) as f:
        download_date = f.readline().split(":")[0]
if snakemake.params.update_db and download_date == today and os.path.exists(snakemake.params.mlst_db):
    log_text = "{}: Database already downloaded today, skipping download.".format(today)
elif snakemake.params.update_db:
    if download_date is None:
        log_text = "{}: no database found. Downloading database".format(today)
    else:
        log_text = "{}: found database updated {}, updating.".format(today, download_date)
    subprocess.Popen("pyngoST.py -d -n {} -cc {}".format(snakemake.params.mlst_db, snakemake.params.ngstar_cc), shell=True).wait()
    with open(download_date_file, 'w') as f:
        f.write("{}: database downloaded.".format(today))
elif download_date is None or not os.path.exists(snakemake.params.mlst_db):
    raise ValueError("Database file missing, please point to pyngost database or set update_db to True.")
else:
    log_text = "Database update not requested: using database downloaded on {}.".format(download_date)

# Create rplF pymlst database
subprocess.Popen("claMLST create --force {mlst_db}/pymlst/pymlst_rplf {mlst_db}/rplf/rplf.txt {mlst_db}/rplf/rplF.tfa".format(
    mlst_db=snakemake.params.mlst_db), shell=True).wait()

# create MLST pymlst database
with open("{}/MLST_profiles.tab".format(snakemake.params.mlst_db)) as f, open("{}/pymlst/MLST_profiles_pymlst.tab".format(snakemake.params.mlst_db), 'w') as o:
    for line in f:
        o.write("\t".join(line.split("\t")[:8]) + "\n")
subprocess.Popen("claMLST create --force {mlst_db}/pymlst/pymlst_mlst {mlst_db}/pymlst/MLST_profiles_pymlst.tab {mlst_db}/abcZ.fas "
    "{mlst_db}/adk.fas {mlst_db}/aroE.fas {mlst_db}/fumC.fas {mlst_db}/gdh.fas {mlst_db}/pdhC.fas {mlst_db}/pgm.fas".format(
    mlst_db=snakemake.params.mlst_db), shell=True).wait()

# create NG-STAR database
with open("{}/NGSTAR_profiles.tab".format(snakemake.params.mlst_db)) as f, open("{}/pymlst/NGSTAR_profiles_pymlst.tab".format(snakemake.params.mlst_db), 'w') as o:
    o.write("ST")
    for i in f.readline().split("\t")[1:8]:
        o.write("\t{}_pymlst".format(i))
    o.write("\n")
    for line in f:
        ST = line.split()[0]
        penA = "{:.3f}".format(float(line.split()[1])).replace('.', '')
        o.write("\t".join([ST, penA] + line.split()[2:-1]) + "\n")
with open("{}/penA.fas".format(snakemake.params.mlst_db)) as f, open("{}/pymlst/penA_pymlst.fas".format(snakemake.params.mlst_db), 'w') as o:
    for line in f:
        if line.startswith(">"):
            o.write(">penA_{}\n".format("{:.3f}".format(float(line.rstrip().split("_")[1])).replace('.', '')))
        else:
            o.write(line)
for i in ["mtrR", "porB", "ponA", "gyrA", "parC", "23S"]:
    with open("{}/{}.fas".format(snakemake.params.mlst_db, i)) as f, open("{}/pymlst/{}_pymlst.fas".format(snakemake.params.mlst_db, i), 'w') as o:
        for line in f:
            if line.startswith(">"):
                o.write(line.split('.')[0] + "\n")
            else:
                o.write(line)

subprocess.Popen("claMLST create --force {mlst_db}/pymlst/pymlst_ngstar {mlst_db}/pymlst/NGSTAR_profiles_pymlst.tab {mlst_db}/pymlst/penA_pymlst.fas "
    "{mlst_db}/pymlst/mtrR_pymlst.fas {mlst_db}/pymlst/porB_pymlst.fas {mlst_db}/pymlst/ponA_pymlst.fas {mlst_db}/pymlst/gyrA_pymlst.fas {mlst_db}/pymlst/parC_pymlst.fas {mlst_db}/pymlst/23S_pymlst.fas".format(
    mlst_db=snakemake.params.mlst_db), shell=True).wait()


# create NG-MAST pymlst database
with open("{}/NGMAST_profiles.tab".format(snakemake.params.mlst_db)) as f, open("{}/pymlst/NGMAST_profiles_pymlst.tab".format(snakemake.params.mlst_db), 'w') as o:
    o.write("ST\tPOR\tTBPB\n")
    f.readline()
    for line in f:
        o.write(line)
subprocess.Popen("claMLST create --force {mlst_db}/pymlst/pymlst_ngmast {mlst_db}/pymlst/NGMAST_profiles_pymlst.tab {mlst_db}/POR.fas "
                 "{mlst_db}/TBPB.fas".format(mlst_db=snakemake.params.mlst_db), shell=True).wait()



with open(snakemake.output.updated_db, 'w') as o, open(snakemake.input.rplf) as rplf_log:
    o.write("rplf updated: " + rplf_log.readline())
    o.write(log_text)