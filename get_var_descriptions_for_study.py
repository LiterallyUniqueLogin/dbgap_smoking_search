#!/usr/bin/env python3

import argparse
import json

import Bio.Entrez

# docs for using Bio.Entrez
# Bio.Entrez docs https://biopython.org/docs/dev/Tutorial/chapter_entrez.html
# Example esearch query in the web interface https://www.ncbi.nlm.nih.gov/gap/?term=2[s_discriminator]%20AND%20phs000280&report=SVariables
#   we want the variable descriptions and IDs from this programmatically
#   maybe also N and variable short names
# Entrez utility descriptions https://www.ncbi.nlm.nih.gov/books/NBK25499/

# was previously using the dbGaP advanced search
# e.g. https://www.ncbi.nlm.nih.gov/gap/advanced_search/?TERM=phs001514
# but "save results" for too many variables causes it to crash
# so we're finding another way

parser = argparse.ArgumentParser()
parser.add_argument('email')
parser.add_argument('phs')
args = parser.parse_args()

Bio.Entrez.email = args.email

# 2[s_discriminator] seems to mean variables
# see here: https://www.ncbi.nlm.nih.gov/gap/?term=2[s_discriminator]%20AND%20phs000280&report=SVariables
search_stream = Bio.Entrez.esearch(db='gap', term=f'2[s_discriminator] AND {args.phs}', usehistory='y')
search_record = Bio.Entrez.read(search_stream)
query_key = search_record['QueryKey']
web_env = search_record['WebEnv']
n_results = int(search_record['Count'])
#n_results = 1000 

#ids = search_record['IdList']
#
#post_stream = Bio.Entrez.epost('gap', id=','.join(ids))
#post_record = Bio.Entrez.read(post_stream)
#query_key = post_record['QueryKey']
#web_env = post_record['WebEnv']

batch_size = 500 # can't download more than 500 in json format at once
# and Bio.Entrez.read() fails for esummary stream XML for whatever reason, so have to use json format
# unless we want to parse XML ourselves
phenotypes = {}
for batch_start in range(0, n_results, batch_size):
    summary_stream = Bio.Entrez.esummary(
        db='gap',
        query_key=query_key,
        webenv=web_env,
        retstart=batch_start,
        retmax=batch_size,
        retmode='json'
    )
    summary_json = json.loads(summary_stream.read())['result']
    for uid in summary_json['uids']:
        uid_result = summary_json[uid]
        assert uid_result['d_object_type'] == 'variable'
        uid_result = uid_result['d_variable_results']

        variable_id_split = uid_result['d_variable_id'].split('|')
        assert len(variable_id_split) == 2
        assert variable_id_split[1][:3] == 'phv'
        variable_id = variable_id_split[1]

        name = uid_result['d_variable_name']
        description = uid_result['d_variable_description'] 
        if any(word in description.lower() for word in ['nicotin', 'smoke', 'smoki', 'cigar', 'vape', 'vapi', 'tobacco', 'hookah']):
            phenotypes[variable_id] = f'{name}\t{variable_id}\t{description}'

for variable_id in sorted(phenotypes):
    print(phenotypes[variable_id])
