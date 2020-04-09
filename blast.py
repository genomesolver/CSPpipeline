import requests
import re
import time
import csv

# ## Initialize query parameters
# ## Check the parameter list at https://ncbi.github.io/blast-cloud/dev/api.html
Forward_Filter = 'ENTREZ_QUERY=txid2[ORGN]'  # Limit to Bacteria
Reverse_Filter = 'ENTREZ_QUERY=txid10239[ORGN]'  # Limit to Viruses
url_endpoint = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?'
csv_header_row = 'query_id,subject_id,per_identity,alignment_length,mismatches,gap_opens,q_start,q_end,s._start,' \
                 's._end,evalue,bit_score,sequence\n'

# ## ------- Auxiliary Routines and Classes -----------------
# A simple routine that extracts an Attribute from a Context given the surrounding marking strings
# ##
def extract_attribute(data, attribute):
    # print("Looking for the attribute:", attribute)
    for line in data.splitlines():
        if attribute in line:
            #print(line)
            attribute_value = re.sub(attribute, "", line)
            attribute_value = re.sub(r"\s", "", attribute_value)
            return attribute_value


# ## ------- End Auxiliary Routines and Classes -----------------

# ## ------- Auxiliary Routines and Classes -----------------
# A simple routine that checks the status of the given RID
# ##

def check_request_status(requestid):
    print("Checking status of RID:", requestid)
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


def blast(protein, query_filter, output_file_name):
    ###############
    # STEP 1 - Submit the query
    ###############

    url_request = 'QUERY=' + protein + '&DATABASE=nr&PROGRAM=blastp&' + query_filter + '&CMD=Put'
    url_submit = url_endpoint + url_request
    #print(url_submit)
    print (output_file_name)

    # Submit the request to the BLAST site
    Submit_Request = requests.put(url_submit)

    # Save the Submit result for troubleshooting
    f = open("tmp/query-submit.html", "w")
    f.write(Submit_Request.text)
    f.close()

    rid = extract_attribute(Submit_Request.text, "RID = ")
    #print(rid)
    rtoe = extract_attribute(Submit_Request.text, "RTOE = ")
    #print(rtoe)

    ###############
    # STEP 2 - Check the status
    ###############
    # https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=3BZF6X4G016

    while True:
        Status, Hits = check_request_status(rid)
        if Status == 'READY' and Hits == 'yes':
            print("We got some hits yay ! ok will grab the hits file")
            break
        elif Status == 'READY' and Hits == 'no':
            print("There seems to be a problem, we have no hits")
            print("Check via: https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=", rid)
            break
        else:
            print("Not ready, will wait for 60 sec and check again")
            #print("You can manually check the status @:\n https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=",rid)
            time.sleep(60)

    ###############
    # STEP 3 - Download the results
    ###############

    # Download the CSV File
    # https://blast.ncbi.nlm.nih.gov/Blast.cgi?RESULTS_FILE=on&FORMAT_TYPE=CSV&FORMAT_OBJECT=Alignment&DESCRIPTIONS=10&ALIGNMENT_VIEW=Tabular&CMD=Get&RID=3BV6047Z014
    url_request = 'RESULTS_FILE=on&FORMAT_TYPE=CSV&FORMAT_OBJECT=Alignment&DESCRIPTIONS=11&ALIGNMENT_VIEW=Tabular&CMD=Get' \
                  '&RID=' + rid
    Submit_DownloadRequest = requests.get(url_endpoint + url_request)
    save_file_handle = open(output_file_name, "w")
    save_file_handle.write(csv_header_row)
    save_file_handle.write(Submit_DownloadRequest.text)
    save_file_handle.close()

    # print("Downloaded the file")

    ##
    # STEP 4
    ##
    # Parse the csv file and write the top 10 results

    csv_file = open(output_file_name, "r")
    dict_reader = csv.DictReader(csv_file)
    i = 0
    print("########################################")
    print("Percentage \tSubject ID\tE_Value")
    print("########################################")
    for row in dict_reader:
        if float(row['per_identity']) > 70:
            # Save the top hit
            if i == 1:
                print(row['per_identity'], "\t", row['subject_id'], "\t", row['evalue'], "\t*")
                i += 1
                topHit = row['subject_id']
            else:
                i += 1
                print(row['per_identity'], "\t", row['subject_id'], "\t", row['evalue'])
    print("########################################")
    print("The top hit = ", topHit)
    return topHit


# ## ------- End Auxiliary Routines and Classes -----------------

# ##
# MAIN PROGRAM
# ##

# Input file with a list of accession numbers
input_file = open("input.txt", "r")
for accession_number in input_file:
    print("Working on the Accession Number:", accession_number.rstrip())
    # Forward blast the accession number
    file_name_slug = accession_number.replace(".1", "")
    file_name = "output/" + file_name_slug.rstrip() + "-results.csv"
    topHit = blast(accession_number, Forward_Filter, file_name)
    # Reverse blast the top hit
    if topHit:
        print("We have a top hit, running a reverse blast on the Accession Number:", topHit)
        reverse_file_name_slug = topHit.replace(".1", "")
        file_name = "output/" + file_name_slug.rstrip() + "-reverse-" + reverse_file_name_slug.rstrip() + "-results.csv"
        blast(topHit, Reverse_Filter, file_name)
input_file.close()


## Troubleshooting (Replace RID with the actual request ID)
## Status:
## https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=RID
## Results:
## https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=RID