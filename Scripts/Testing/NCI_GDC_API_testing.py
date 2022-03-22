import requests
import json
import re

files_endpt = "https://api.gdc.cancer.gov/files"

filters = {
    "op": "and",
    "content":[
        {
        "op": "in",
        "content":{
            "field": "cases.project.primary_site",
            "value": ["Lung"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "cases.demographic.race",
            "value": ["white"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "cases.demographic.gender",
            "value": ["female"]
            }
        },
        {
        "op": "in",
        "content":{
            "field": "files.analysis.workflow_type",
            "value": ["HTSeq - FPKM"]
            }
        }
    ]
}

# Here a GET is used, so the filter parameters should be passed as a JSON string.

params = {
    "filters": json.dumps(filters),
    "fields": "file_id",
    "format": "JSON",
    "size": "1000"
    }

response = requests.get(files_endpt, params = params)

file_uuid_list = []

# This step populates the download list with the file_ids from the previous query
for file_entry in json.loads(response.content.decode("utf-8"))["data"]["hits"]:
    file_uuid_list.append(file_entry["file_id"])

data_endpt = "https://api.gdc.cancer.gov/data"

params = {"ids": file_uuid_list}

response = requests.post(data_endpt, data = json.dumps(params), headers = {"Content-Type": "application/json"})

response_head_cd = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head_cd)[0]

with open(file_name, "wb") as output_file:
    output_file.write(response.content)


# Script for testing to get data from National Cancer Institure using the GDCs API
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

