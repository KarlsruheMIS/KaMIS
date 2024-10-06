#!/bin/bash

run="mmwis_soa_full_graph"
program="mmwis"
run_parallel="run_parallel_${run}.txt"

path_to_program="./deploy/${program}"
path_to_instance="/home/ernestineg/test_instances"

create_run() {
    id=$1
    id_instructions=$2
    result_file=$3
    error_file=$4
    address2=${address}/${id}
    mkdir -p ${address2}
    
    instruction="${path_to_program} ${id_instructions}"
    echo " ${instruction} > ${address2}/result.txt || { echo '${id} ${address2}' >> ${error_file}; echo 'error,' >> ${result_file}; } && { res='${id}'\`cat ${address2}/result.txt | grep -E 'Time found:|Weight:' | awk '{printf \",\" \$NF}'\` && echo \${res}  >> ${result_file}; }; " >> ${run_parallel}
}
rm -f ${run_parallel}


    #########################################################
    ############### Experiments for parameter tuning ########
    #########################################################
    instance="chils_set"
    address="results/${run}"
    database="${path_to_instance}/${instance}"
    results="results/results_${run}.csv"
    errors="results/error_${run}.txt"
    echo "graph,weight,time" > ${results}
    mkdir -p $address
    for graph in `ls ${database} | grep "graph" | awk -F"." '{print $1}'`; do
        overall=" --config=mmwiss ${database}/${graph}.graph --time_limit 3600"
        create_run "${graph}" "${overall}" "${results}" "${errors}"
    done
