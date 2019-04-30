#!/bin/bash
#PBS -lwalltime=---job_dur---
#PBS -lnodes=1

cd $HOME/projects/retino_HCP/
python ---fit_file--- ---fit_model--- ---subject--- ---start_idx--- ---end_idx--- ---data_file--- ---base_dir---