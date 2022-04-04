# Imports
import sys

import requests
import json
import re
import os, glob
import tarfile
import gzip
from pathlib import Path
import shutil


def main():
    project_dir = Path.cwd()

    file_name = getCountData([''], project_dir)
    extractFiles(file_name, project_dir, True)
    return None


def getCountData(user_input, project_dir):
    """
    #TODO Add documentation
    :param user_input:
    :param project_dir:
    :return:
    """
    # TODO: Add user paramenter wildcards - see Git Joinks.
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
                    # "value": ["HTSeq - Counts"]  # This seems to be removed from GDC
                    "value": ["STAR - Counts"]
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

    # Check if there's items in the supplemented list.
    file_entries = json.loads(response.content.decode("utf-8"))["data"]["hits"]
    if file_entries:
        print('Sucessfully found files, interpreting...')
    else:
        print('No data was found with the given search criteria!')
        print('Exiting now...')
        sys.exit(0)

    # This step populates the download list with the file_ids from the previous query
    counter = 0
    for file_entry in file_entries:
        file_uuid_list.append(file_entry["file_id"])
        counter += 1
        # if counter >= 10: # Debugging
        #     break
    print('Total number of entries found:\t', counter)
    print('Fetching data files from entries...')

    params = {"ids": file_uuid_list}
    response = requests.post(data_endpt, data=json.dumps(params), headers={"Content-Type": "application/json"})
    response_head_cd = response.headers["Content-Disposition"]
    file_name = re.findall("filename=(.+)", response_head_cd)[0]

    print('Writing')
    with open(Path(project_dir / 'Output' / 'Counts' / file_name), "wb") as output_file:
        output_file.write(response.content)
    return file_name


def extractFiles(file_name, project_dir, remove_zip=True):
    print('Exracting main tar file')
    my_tar = tarfile.open(Path(project_dir / 'Output' / 'Counts' / file_name))
    out_file = file_name.rsplit('.tar.gz', 1)[0]  # Can't use stem since there's two extions behind it .tar.gz
    my_tar.extractall(Path(project_dir / 'Output' / 'Counts' / out_file))  # specify which folder to extract to
    my_tar.close()

    if remove_zip:
        os.remove(Path(project_dir / 'Output' / 'Counts' / file_name))

    name = ''
    print('Unpacking...')
    for root, dirs, files in os.walk(Path(project_dir / 'Output' / 'Counts' / out_file)):
        output_dir = Path(root).parent
        for name in files:

            if name.endswith(".tsv"):
                name = Path(name)
                shutil.move(Path(root / name), Path(output_dir / name))
                os.removedirs(Path(root))

    print('Successfully unpacked '+str(name)+ '& removed zipped files!')
    print('Please note:\tFile name consists of projectID _ fileID of GDC database.')
    return None

def decompress(infile, tofile):
    with open(infile, 'rb') as inf, open(tofile, 'w', encoding='utf8') as tof:
        decom_str = gzip.decompress(inf.read()).decode('utf-8')
        tof.write(decom_str)


# Note: This might need to be changed when adding the scripts to a pipeline
if __name__ == "__main__":
    main()