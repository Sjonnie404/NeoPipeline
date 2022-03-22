# Script for testing to get data from National Cancer Institure using the GDCs API

### New #####
import requests
import json
import re
import os

project_path = os.getcwd()
pp = project_path
pp = pp.rsplit('\\', 1)[0]


files_endpt = "https://api.gdc.cancer.gov/files"

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
            "value": ["female"]
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


data_endpt = "https://api.gdc.cancer.gov/data"

params = {"ids": file_uuid_list}

response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type": "application/json"})

response_head_cd = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head_cd)[0]

print('Writing')
with open(pp+'\\Output\\Counts\\'+file_name, "wb") as output_file:
    output_file.write(response.content)


# PLEASE NOTE: DEPRECATED

# import requests
# import json
# import re
#
# file_id = "b658d635-258a-4f6f-8377-767a43771fe4"
#
# data_endpt = "https://api.gdc.cancer.gov/data/{}".format(file_id)
#
# response = requests.get(data_endpt, headers = {"Content-Type": "application/json"})
#
# # The file name can be found in the header within the Content-Disposition key.
# response_head_cd = response.headers["Content-Disposition"]
#
# file_name = re.findall("filename=(.+)", response_head_cd)[0]
#
# with open(file_name, "wb") as output_file:
#     output_file.write(response.content)
