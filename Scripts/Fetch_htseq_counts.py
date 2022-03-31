# Script for testing to get data from National Cancer Institure using the GDCs API

### New #####
import shutil

import requests
import json
import re
import os, glob
import tarfile
import gzip

project_path = os.getcwd()
pp = project_path
pp = pp.rsplit('\\', 1)[0]


def main():
    getCountData()
    extractFiles()
    return None


def getCountData(user_input=''):
    files_endpt = "https://api.gdc.cancer.gov/files"
    data_endpt = "https://api.gdc.cancer.gov/data"
    filters = {
        "op": "and",
        "content":[
            {
            "op": "in",
            "content":{
                "field": "cases.project.primary_site",
                "value": ["Skin"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": ["TCGA-SKCM"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.analysis.workflow_type",
                    "value": ["HTSeq - Counts"]
                }
            },
            {
            "op": "in",
            "content":{
                "field": "files.data_category",
                "value": ["transcriptome profiling"]
                }
            },
            {
            "op": "in",
            "content":{
                "field": "files.data_type",
                "value": ["Gene Expression Quantification"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.experimental_strategy",
                    "value": ["RNA-Seq"]
                }
            }
        ]
    }

    # Here a GET is used, so the filter parameters should be passed as a JSON string.
    params = {
        "filters": json.dumps(filters),
        "fields": "file_id",
        "format": "JSON",
        "size": '9999999'  # For some reason there's no all function like 0 or -1, and larger numbers break
        }

    response = requests.get(files_endpt, params=params)
    file_uuid_list = []
    print('Started populating download list...')

    # This step populates the download list with the file_ids from the previous query
    index = 0
    for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
        file_uuid_list.append(file_entry["file_id"])
        index += 1
        print('added entry: '+str(index))

    params = {"ids": file_uuid_list}
    response = requests.post(data_endpt, data=json.dumps(params), headers={"Content-Type": "application/json"})
    response_head_cd = response.headers["Content-Disposition"]
    file_name = re.findall("filename=(.+)", response_head_cd)[0]

    exit()
    print('Writing')
    with open(pp+'\\Output\\Counts\\'+file_name, "wb") as output_file:
        output_file.write(response.content)
    return None


# TODO: Add  variable for wd parth and output folder path
def extractFiles(remove_zip=True):
    file_name=''

    print('Exracting main tar file')
    my_tar = tarfile.open(pp+'\\Output\\Counts\\'+file_name)
    out_file = file_name.rsplit('.tar.gz', 1)[0]
    my_tar.extractall(pp+'\\Output\\Counts\\'+out_file+'\\')  # specify which folder to extract to
    my_tar.close()
    os.remove(pp+'\\Output\\Counts\\'+file_name)

    print('Unpacking...')
    for root, dirs, files in os.walk(pp+'\\Output\\Counts\\'+out_file+'\\'):
        for name in files:
            if name.endswith(".counts.gz"):
                dir_name = root.rsplit('\\', 2)
                subdir_name = dir_name[2]
                dir_name = dir_name[1]

                decompress(root+'\\'+name, pp+'\\Output\\Counts\\'+dir_name+'\\'+subdir_name+'_'+name.rsplit('.counts.gz',1)[0])

                os.remove(root+'\\'+name)
                os.removedirs(root)
    print('Successfully unpacked & removed zipped files!')
    print('Please note:\tFile name consists of projectID _ fileID of GDC database.')
    return None


def decompress(infile, tofile):
    with open(infile, 'rb') as inf, open(tofile, 'w', encoding='utf8') as tof:
        decom_str = gzip.decompress(inf.read()).decode('utf-8')
        tof.write(decom_str)


# Note: This might need to be changed when adding the scripts to a pipeline
if __name__ == "__main__":
    main()