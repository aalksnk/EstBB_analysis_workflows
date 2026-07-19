#!/bin/sh

original_path=[Original location of all .parquet files]
pipeline_input_path=[Input folder for pipeline, where symlinks are made]

mkdir -p ${pipeline_input_path}

paths=$(ls ${original_path})
#If it is needed to input only subset of files (i.e. 100):
#paths=$(ls ${original_path} | head -n 100)

for file in ${paths}
do 
ln -sf ${original_path}/${file} ${pipeline_input_path}/${file}
done