#!/bin/bash
echo "target_name,target_accession,target_length,query_name,query_accession,query_length,evalue,score,bias,domain,ndomains,domain_c_evalue,domain_i_evalue,domain_score,domain_bias,hmm_from,hmm_to,ali_from,ali_to,env_from,env_to,acc,description"
grep -v "^#" $1 | perl -pe 's/ +/,/g'
