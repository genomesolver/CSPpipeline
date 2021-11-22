import pandas as pd
import requests
import re
import time

start_time = time.time() # Change: begin measuring runtime

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

# A small function to convert percentages from strings to floats
def str_to_float(string):
    return float(string.strip('%'))

# ## ------- End Auxiliary Routines and Classes -----------------

# ## ------- Auxiliary Routines and Classes -----------------
# A simple routine that checks the status of the given RID
# ##

def check_request_status(requestid):
    # print("Checking status of RID:", requestid)
    # Submit the request to the BLAST site
    submit_request = requests.put('https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=' + requestid) 

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


def blast(protein, query_filter):
    ###############
    # STEP 1 - Submit the query
    ###############

    # Submit the request to the BLAST site
    Submit_Request = requests.put('https://blast.ncbi.nlm.nih.gov/Blast.cgi?QUERY=' + protein + '&DATABASE=nr&PROGRAM=blastp&' + query_filter + '&CMD=Put') 

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
    status, hits = '', '' # Initialize status and hits variables
    while status != 'READY': 
        status, hits = check_request_status(rid)
        if status != 'READY': # If updated status still not ready, wait one minute
            print("Will wait and check in 60 seconds")  # debug msg
            time.sleep(60) 
    if hits == 'yes': # At this point status must be ready, check if there are hits 
        print("Search successful! Selecting top hit from results table")
    elif hits == 'no': # Print error message and return None if there were no hits
        print("Error: The search did not return any hits. Please try again with adjusted parameters or different accession numbers. See https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=" + rid + " for more details.")
        return None 

    ###############
    # STEP 3
    ###############

    # Use pandas web scraping to turn results table into dataframe
    # Process dataframe for the top hit
    url = 'https://blast.ncbi.nlm.nih.gov/Blast.cgi?ADV_VIEW=on&ALIGNMENTS=100&ALIGNMENT_VIEW=Pairwise&RESULTS_FILE=on&DYNAMIC_FORMAT=on&FORMAT_TYPE=HTML&FORMAT_OBJECT=Alignment&CMD=Get&RID=' + rid
    hit_table = pd.read_html(url)[4]
    hit_table['Query Cover'] = hit_table['Query Cover'].apply(str_to_float) # Change query cover % and % identity to floats
    hit_table['Per. Ident'] = hit_table['Per. Ident'].apply(str_to_float)
    hit_table = hit_table.loc[hit_table['Query Cover']>70].iloc[:1] # Select the top hit whose query cover % is greater than 70
    hit_table = hit_table.filter(items=['Scientific Name','Query Cover','E value','Per. Ident','Accession']) # Select the relevant columns
    topHit = hit_table.iloc[0]
    print(hit_table.to_string()) # Print results dataframe of top hit to screen
    print("###################################################")
    print("##### TOP HIT = " + topHit['Accession'], topHit['Scientific Name'])
    print("###################################################")
    return topHit['Accession'] # Return top hit accession number


# ## ------- End Auxiliary Routines and Classes -----------------

# #################
#   MAIN PROGRAM
# #################

# Input file with a list of accession numbers
input_file = open("input.txt", "r")
save_csv_file_name = 'output/all-results.csv' 
queries = []
hits = []
for accession_number in input_file:
    queries.append(accession_number)
    print("Working on the Accession Number:", accession_number.rstrip())
    # Forward blast the accession number
    result = blast(accession_number, 'ENTREZ_QUERY=txid2[ORGN]')
    # Reverse blast the top hit's accession number
    hits.append(result)
    if result is not None:
        print("We have a top hit, running a reverse blast on the Accession Number:", result)
        queries.append(result)
        rev_result = blast(result, 'ENTREZ_QUERY=txid10239[ORGN]')
        if rev_result is not None:
            hits.append(rev_result)
        else:
            hits.append('No result')
    else:
        queries.append('No result')
        hits.append('No result')
input_file.close()
if rev_result == accession_number:
    print('The search has returned a positive result indicating horizontal gene transfer.\nForward Blast: {} -> {} Reverse Blast: {} -> {}'.format(accession_number, result[0], result[0], rev_result[0]))
data_dict = {'Query ID':queries, 'Search Results':hits}
results_df = pd.DataFrame(data_dict)
results_df.to_csv(path_or_buf=save_csv_file_name)
runtime = (time.time()-start_time) # Change: calculate and store total program runtime
print("Runtime: {} seconds".format(runtime)) # Change: print runtime to screen

## Troubleshooting (Replace RID with the actual request ID)
## Status:
## https://blast.ncbi.nlm.nih.gov/Blast.cgi?CMD=Get&FORMAT_OBJECT=SearchInfo&RID=RID
## Results:
## https://blast.ncbi.nlm.nih.gov/blast/Blast.cgi?CMD=Get&FORMAT_TYPE=Text&RID=RID