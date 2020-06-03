#!/bin/bash

# required parameters
handle='test'
data_filename='data1.txt'
box_size=1024
corr_type='CF_monopole'
dim1_min=0
dim1_max=120
dim1_nbin=60

# only used for cross-correlation
data_filename_2=''

# only used for r-mu or s-pi binning
dim2_min=0
dim2_max=0
dim2_nbin=0

python $HOME/code/CONTRAST/contrast.py \
--handle "$handle" \
--data_filename "$data_filename" \
--data_filename_2 "$data_filename_2" \
--box_size "$box_size" \
--corr_type "$corr_type" \
--dim1_min "$dim1_min" \
--dim1_max "$dim1_max" \
--dim1_nbin "$dim1_nbin" \
--dim2_min "$dim2_min" \
--dim2_max "$dim2_max" \
--dim2_nbin "$dim2_nbin"


