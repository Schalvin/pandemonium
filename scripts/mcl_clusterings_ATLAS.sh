#!/bin/bash

ADJACENCY_TABLE="/home/chili/chalvin/projet_long/ATLAS/ATLAS_SWORD2_filt_adjacency.tsv"

mcl $ADJACENCY_TABLE --abc -I 1.35
mcl $ADJACENCY_TABLE --abc -I 1.5
mcl $ADJACENCY_TABLE --abc -I 1.7
mcl $ADJACENCY_TABLE --abc -I 2
mcl $ADJACENCY_TABLE --abc -I 3
mcl $ADJACENCY_TABLE --abc -I 4.5