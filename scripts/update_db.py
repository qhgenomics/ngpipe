import os
from urllib.request import urlopen
import json

db = snakemake.params.db
if db == "ngstar":
    scheme = "67"
elif db == "ngmast":
    scheme = "71"
elif db == "mlst":
    scheme = "1"

params_mlst_dir = snakemake.params.mlst_dir

dbdir = os.path.join(params_mlst_dir, 'db', 'pubmlst', db)
version_log = os.path.join(dbdir, "version.log")
if not os.path.exists(dbdir):
    os.makedirs(dbdir)
if os.path.exists(version_log):
    with open(version_log) as f:
        version = f.readline().rstrip()
else:
    version = None


with urlopen("https://rest.pubmlst.org/db/pubmlst_neisseria_seqdef/schemes/" + scheme) as response:
    response_content = response.read().decode('utf-8')
scheme_info = json.loads(response_content)


online_version = scheme_info["last_added"]

if version != online_version:
    profiles = scheme_info["profiles_csv"]
    response = urlopen(profiles)
    data = response.read()
    profiles_local = os.path.join(dbdir, db + '.txt')
    with open(profiles_local, 'wb') as o:
        o.write(data)
    for i in scheme_info['loci']:
        loci = i.split('/')[-1]
        response = urlopen(i)
        loci_info = json.loads(response.read().decode('utf-8'))
        loci_fasta = loci_info["alleles_fasta"]
        response = urlopen(loci_fasta)
        data = response.read()
        loci_local = os.path.join(dbdir, loci + '.tfa')
        loci_local = loci_local.replace("'", "")
        loci_local = loci_local.replace("NG-MAST_", "")
        with open(loci_local, 'wb') as o:
            o.write(data)
        if db == "ngstar":
            loci_alleles = loci_info["alleles"]
            loci_alleles += "?return_all=1"
            response = urlopen(loci_alleles)
            json_alleles = json.loads(response.read())
            with open(loci_local + '.comments', 'w') as o:
                for j in json_alleles["alleles"]:
                    response = urlopen(j)
                    json_allele = json.loads(response.read())
                    allele_id = json_allele["allele_id"]
                    if "comments" in json_allele:
                        comments = json_allele["comments"]
                    else:
                        comments = "None"
                    o.write("{}\t{}\n".format(allele_id, comments))


    with open(version_log, 'w') as o:
        o.write(online_version + "\n")
    with open(snakemake.output.log, 'w') as o:
        o.write("Database updated to version uploaded on " + online_version + "\n")
else:
    with open(snakemake.output.log, 'w') as o:
        o.write("Database already using version uploaded on " + online_version + "\n")