#!/bin/bash

dockername="$2"

if [[ -z "$dockername" ]]; then
    dockername="bio_condensates"
fi

docker stop $dockername && docker remove $dockername
docker build -f dockerfile . -t $dockername:latest

docker run \
    --name=$dockername \
    -v /home/max/projects/bio_condensates/output:/bio_condensates/output \
    $dockername:latest \
    $1
