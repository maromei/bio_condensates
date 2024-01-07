#!/bin/bash

usage() {
    echo "
    Usage: $0
        [-n sim_name]
        [-d dockername=bio_condensates]
        [-o output_dir=/home/max/projects/bio_condensates/output]
        [-i init_file=init/{sim_name}.dat]
    ";
}

while getopts ":n:d:o:i:" opt; do
    case "${opt}" in
        n) sim_name_to_run="${OPTARG}"
        ;;
        d) dockername="${OPTARG}"
        ;;
        o) output_dir="${OPTARG}"
        ;;
        i) init_file="${OPTARG}"
        ;;
        \?) usage; exit -1;
        ;;
    esac
done

if [[ -z "$sim_name_to_run" ]]; then
    echo "You need to choose a simulation name to run"
    usage
    exit -1
fi

if [[ -z "$dockername" ]]; then
    dockername="bio_condensates"
fi

if [[ -z "$output_dir" ]]; then
    output_dir="/home/max/projects/bio_condensates/output"
fi

if [[ -z "$init_file" ]]; then
    init_file="init/$sim_name_to_run.dat"
fi

docker stop $dockername && docker remove $dockername
docker build -f dockerfile . -t $dockername:latest

docker run \
    --name=$dockername \
    -v $output_dir:/bio_condensates/output \
    $dockername:latest \
    $sim_name_to_run
    $init_file
