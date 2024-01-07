#!/bin/bash

sim_name_to_run="$1"
dockername="$2"
output_dir="$3"

if [[ -z "$dockername" ]]; then
    dockername="bio_condensates"
fi

if [[ -z "$output_dir" ]]; then
    output_dir="/home/max/projects/bio_condensates/output"
fi

docker stop $dockername && docker remove $dockername
docker build -f dockerfile . -t $dockername:latest

docker run \
    --name=$dockername \
    -v $output_dir:/bio_condensates/output \
    $dockername:latest \
    $sim_name_to_run
