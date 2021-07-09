import requests
import re
import time
import json

# ## Initialize query parameters
# ## Check the parameter list at https://ncbi.github.io/blast-cloud/dev/api.html
Forward_Filter = 'ENTREZ_QUERY=txid2[ORGN]'  # Limit to Bacteria
Reverse_Filter = 'ENTREZ_QUERY=txid10239[ORGN]'  # Limit to Viruses
url_endpoint = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?'
output_csv_header_row = 'query_id,scientific_name,query_cover_per,evalue,per_identity,accession_id\n'


# ## ------- Auxiliary Routines and Classes -----------------
# A simple routine that extracts an Attribute from a Context given the surrounding marking strings
# ##
def extract_attribute(data, attribute):
    # print("Looking for the attribute:", attribute)
    for line in data.splitlines():
        if attribute in line:
            # print(line)
            attribute_value = re.sub(attribute, "", line)
            attribute_value = re.sub(r"\s", "", attribute_value)
            return attribute_value


# ## ------- End Auxiliary Routines and Classes -----------------

# ## ------- Auxiliary Routines and Classes -----------------
# A simple routine that checks the status of the given RID
# ##

def check_request_status(requestid):
    # print("Checking status of RID:", requestid)
    url_request = 'CMD=Get&FORMAT_OBJECT=SearchInfo&RID=' + requestid
    url_submit = url_endpoint + url_request
    # Submit the request to the BLAST site
    submit_request = requests.put(url_submit)

    # Save the Submit result for troubleshooting
    file_handle = open("tmp/query-status.html", "w")
    file_handle.write(submit_request.text)
    file_handle.close()

    query_status = extract_attribute(submit_request.text, "Status=")
    # print(query_status)
    query_hits = extract_attribute(submit_request.text, "ThereAreHits=")
    # print(query_hits)
    return query_status, query_hits


# ## ------- End Auxiliary Routines and Classes -----------------

# ## ------- Auxiliary Routines and Classes -----------------
# Main subroutine to blast a protein forward / reverse
# ##


def blast(protein, query_filter, save_json_file_name, save_csv_file_name):
    ###############
    # STEP 1 - Submit the query
    ###############

    url_request = 'QUERY=' + protein + '&DATABASE=nr&PROGRAM=blastp&' + query_filter + '&CMD=Put'
    url_submit = url_endpoint + url_request

    # Submit the request to the BLAST site
    Submit_Request = requests.put(url_submit)

    # Save the Submit result for troubleshooting
    f = open("tmp/query-submit.html", "w")
    f.write(Submit_Request.text)
    f.close()

    # Extract the request id that will be used for next steps
    rid = extract_attribute(Submit_Request.text, "RID = ")

    ###############
    # STEP 2
    ###############
    print("Check the status via web-browser:")  # debug msg
    print("https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + rid)  # debug msg
    while True:
        Status, Hits = check_request_status(rid)
        if Status == 'READY' and Hits == 'yes':
            print("We got some hits, will grab the json hits file")
            break
        elif Status == 'READY' and Hits == 'no':
            print("Error: We have no hits")  # debug msg
            print(
                "Check: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + rid)
            break
        else:
            print("Will wait and check in 60 seconds")  # debug msg
            time.sleep(60)

    ###############
    # STEP 3
    ###############
    # Download the JSON File to the tmp directory
    # https://blast.ncbi.nlm.nih.gov/Blast.cgi?RESULTS_FILE=on&FORMAT_TYPE=JSON2_S&FORMAT_OBJECT=Alignment&CMD=Get&RID=3BV6047Z014

    url_request = 'RESULTS_FILE=on&FORMAT_TYPE=JSON2_S&FORMAT_OBJECT=Alignment&CMD=Get&RID=' + rid
    Submit_JSONRequest = requests.get(url_endpoint + url_request)
    save_json_file_handle = open(save_json_file_name, "w")
    save_json_file_handle.write(Submit_JSONRequest.text)
    save_json_file_handle.close()

    print("Downloaded the JSON hits file")

    ###############
    # STEP 4
    ###############
    # Parse the JSON file and write the top 10 results to an CSV file in the output folder

    output_csv_file = open(save_csv_file_name, "w")
    # Add the CSV Header row fields
    output_csv_file.write(output_csv_header_row)

    with open(save_json_file_name, "r") as read_file:
        data = json.load(read_file)

    # Extract the query details
    query_id = data['BlastOutput2'][0]['report']['results']['search']['query_id']
    query_len = data['BlastOutput2'][0]['report']['results']['search']['query_len']

    print("###################################################")
    print("Query Cover %\tE_Value\tAccession Id\tSubject Name")
    print("###################################################")

    # We need only the top 10 hits
    hit_count = 0
    for hit in data['BlastOutput2'][0]['report']['results']['search']['hits']:
        # Extract the fields we need
        scientific_name = hit['description'][0]['sciname']
        accession_id = hit['description'][0]['accession']
        hsps_align_len = hit['hsps'][0]['align_len']
        hsps_identity = hit['hsps'][0]['identity']
        hsps_query_from = hit['hsps'][0]['query_from']
        hsps_query_to = hit['hsps'][0]['query_to']
        evalue = hit['hsps'][0]['evalue']

        # How to calculate Query Cover %
        # https://codereview.stackexchange.com/questions/39879/calculate-query-coverage-from-blast-output)
        query_cover_per = ((hsps_query_to - hsps_query_from) / query_len) * 100

        per_identity = (hsps_identity / hsps_align_len) * 100
        if (query_cover_per > 70) and (hit_count < 10):
            if hit_count == 0:
                topHit = accession_id
                topHit_scientific_name = scientific_name
            hit_count += 1
            print(round(query_cover_per, 2), evalue, accession_id, scientific_name, sep='\t')
            row = ','.join(
                (query_id, scientific_name, str(query_cover_per), str(evalue), str(per_identity), accession_id)) + '\n'
            output_csv_file.write(row)
    print("###################################################")
    print("##### TOP HIT = " + topHit, topHit_scientific_name)
    print("###################################################")
    output_csv_file.close()
    return topHit


# ## ------- End Auxiliary Routines and Classes -----------------

# #################
#   MAIN PROGRAM
# #################

# Input file with a list of accession numbers
input_file = open("input.txt", "r")
for accession_number in input_file:
    print("Working on the Accession Number:", accession_number.rstrip())
    # Forward blast the accession number
    file_name_slug = (accession_number.replace(".1", ""))
    save_json_file_name = "tmp/" + file_name_slug.rstrip() + "-results.json"
    save_csv_file_name = "output/" + file_name_slug.rstrip() + "-results.csv"
    topHit = blast(accession_number, Forward_Filter, save_json_file_name, save_csv_file_name)
    # Reverse blast the top hit's accession number
    if topHit:
        print("We have a top hit, running a reverse blast on the Accession Number:", topHit)
        reverse_file_name_slug = topHit.replace(".1", "")
        save_json_file_name = "tmp/" + file_name_slug.rstrip() + "-rev-" + reverse_file_name_slug.rstrip() + "-results.json"
        save_csv_file_name = "output/" + file_name_slug.rstrip() + "-rev-" + reverse_file_name_slug.rstrip() + "-results.csv"
        blast(topHit, Reverse_Filter, save_json_file_name, save_csv_file_name)
input_file.close()

## Troubleshooting (Replace RID with the actual request ID)
## Status:
## https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=RID
## Results:
## https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=RID
