########################################################################################
#  Script to fetch count data from the GDC website.
#  Input: user input on the project and data to fetch
#  Output: Extracted folder with star count data files
#  Made by Shane Pullens, Utrecht University - Theoretical bioinformatics.
# Note, only need to connect to server when using NetMHCpan, and for deployment
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


# TODO: Find a way to automatically download the sample sheet, if not possible make a dummy varient since the data isn't used.
# TODO: Refractor the Output folder to the Data folder
# Note: When selecting a project with a single file, the algorithm resturns a single file instead of a zip file.
# Note: This is not yet supported, so please skip single-file, projects.


def main():
    project_dir = Path.cwd()

    user_input = {'primary_site': 'breast',
                  'project_id': 'EXCEPTIONAL_RESPONDERS-ER'}
    file_name = getCountData(user_input, project_dir)
    extractFiles(file_name, project_dir, user_input, True)

    return None


def getCountData(user_input, project_dir):
    """
    #TODO Add documentation
    :param user_input:
    :param project_dir:
    :return:
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
    # filters = json.loads(filters)

    # Here a GET is used, so the filter parameters should be passed as a JSON string.
    params = {
        "filters": filters,
        "fields": "file_id",
        "format": "JSON",
        "size": '99999'  # For some reason there's no all function like 0 or -1, and larger numbers break
        } #99999

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

    # print(response.content)
    print('Writing...')
 #   print(json.dumps(response.content, indent=2))

    with open(Path(project_dir / 'Output' / 'Counts' / file_name), "wb") as output_file:
        output_file.write(response.content)
    return file_name


def extractFiles(file_name, project_dir, input, remove_zip=True):
    """

    :param file_name:
    :param project_dir:
    :param remove_zip:
    :return:
    """
    print('Exracting main tar file')
    my_tar = tarfile.open(Path(project_dir / 'Output' / 'Counts' / file_name))
    out_file = file_name.rsplit('.tar.gz', 1)[0]  # Can't use stem since there's two extions behind it .tar.gz

    input = input.get('primary_site')+'_'+input.get('project_id')+'_'
    new_out_file = out_file.split('gdc_download_')[1]
    new_out_file = input+new_out_file

    my_tar.extractall(Path(project_dir / 'Output' / 'Counts' / new_out_file))  # specify which folder to extract to
    my_tar.close()

    if remove_zip:
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
                except OSError:  # This throws an error because of non-empty 'raw_counts' directory # TODO Fix this
                    continue
    print('Successfully unpacked '+str(name)+' & removed zipped files!')
    print('Please note:\tFile name consists of Tissue_projectName_projectID_fileID of GDC database.')
    return new_out_file


def decompress(infile, tofile): # NOTE Should be deprecated
    with open(infile, 'rb') as inf, open(tofile, 'w', encoding='utf8') as tof:
        decom_str = gzip.decompress(inf.read()).decode('utf-8')
        tof.write(decom_str)


# Note: This might need to be changed when adding the scripts to a pipeline
if __name__ == "__main__":
    main()