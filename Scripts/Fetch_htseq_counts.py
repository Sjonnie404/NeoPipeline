########################################################################################
#  Script to fetch count data from the GDC website.
#  Input: user input on the project and data to fetch
#  Output: Extracted folder with star count data files
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
########################################################################################
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
    return None


def getCountData(user_input, project_dir):
    """
    Fetches the star count data from the GDC data portal based on selected cancer site & project.
    And writes this to the selected folder
    :param user_input: dictionary that hold the selected cancer site & selected project.
    :param project_dir: defined working folder
    :return: filename of the downloaded tar file.
    """
    files_endpt = "https://api.gdc.cancer.gov/files"
    data_endpt = "https://api.gdc.cancer.gov/data"
    filters = '''{
        "op": "and",
        "content":[
            {
            "op": "in",
            "content":{
                "field": "cases.primary_site",
                "value": ["%s"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "cases.project.project_id",
                    "value": ["%s"]
                }
            },
            {
                "op": "in",
                "content": {
                    "field": "files.analysis.workflow_type",
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
    }''' % (user_input.get('primary_site'), user_input.get('project_id'))

    # Here a GET is used, so the filter parameters should be passed as a JSON string.
    params = {
        "filters": filters,
        "fields": "file_id",
        "format": "JSON",
        "size": '99999'  # For some reason there's no all function like 0 or -1, and larger numbers break
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

    print('Total number of entries found:\t', counter)
    print('Fetching data files from entries...')
    params = {"ids": file_uuid_list}
    response = requests.post(data_endpt, data=json.dumps(params), headers={"Content-Type": "application/json"})
    response_head_cd = response.headers["Content-Disposition"]
    file_name = re.findall("filename=(.+)", response_head_cd)[0]

    print('Writing...')
    with open(Path(project_dir / 'Output' / 'Counts' / file_name), "wb") as output_file:
        output_file.write(response.content)
    return file_name


def extractFiles(file_name, project_dir, input, remove_zip=True):
    """
    Extracts the files from the previously downloaded tar file.
    :param file_name: name of the tarfile
    :param project_dir: name of the working directory
    :param input: user input to define filenames
    :param remove_zip: boolean to remove zip files after extraction
    :return:
    """

    print('Extracting main tar file...')
    my_tar = tarfile.open(Path(project_dir / 'Output' / 'Counts' / file_name))
    # Can't use Path's stem since there's two extensions behind it (.tar.gz)
    out_file = file_name.rsplit('.tar.gz', 1)[0]

    input = input.get('primary_site')+'_'+input.get('project_id')+'_'
    # removes 'gdc_download_'
    new_out_file = out_file.split('gdc_download_')[1]
    new_out_file = input+new_out_file

    my_tar.extractall(Path(project_dir / 'Output' / 'Counts' / new_out_file))  # specify which folder to extract to
    my_tar.close()

    if remove_zip: # When True, remove leftover zip files
        os.remove(Path(project_dir / 'Output' / 'Counts' / file_name))

    name = ''
    print('Unpacking...')
    os.mkdir(Path(project_dir / 'Output' / 'Counts' / new_out_file / 'Raw_counts'))
    for root, dirs, files in os.walk(Path(project_dir / 'Output' / 'Counts' / new_out_file)):
        output_dir = Path(root).parent
        for name in files:
            if name.endswith("counts.tsv"):
                name = Path(name)
                shutil.move(Path(root / name), Path(output_dir / 'Raw_counts' / name))
                try:
                    os.removedirs(Path(root))
                except OSError:  # This throws an error because of non-empty 'raw_counts' directory
                    continue
    print('Successfully unpacked '+str(name)+' & removed zipped files!')
    print('Please note:\tFile name consists of Tissue_projectName_projectID_fileID of GDC database.')
    print(f'Filename variable renamed to:\t{new_out_file}')
    return new_out_file


if __name__ == "__main__":
    main()