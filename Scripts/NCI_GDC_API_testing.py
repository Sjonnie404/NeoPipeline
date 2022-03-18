# Script for testing to get data from National Cancer Institure using the GDCs API

# PLEASE NOTE: DEPRECATED

import requests
import json
import re

file_id = "b658d635-258a-4f6f-8377-767a43771fe4"

data_endpt = "https://api.gdc.cancer.gov/data/{}".format(file_id)

response = requests.get(data_endpt, headers = {"Content-Type": "application/json"})

# The file name can be found in the header within the Content-Disposition key.
response_head_cd = response.headers["Content-Disposition"]

file_name = re.findall("filename=(.+)", response_head_cd)[0]

with open(file_name, "wb") as output_file:
    output_file.write(response.content)
