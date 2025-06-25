#!/bin/bash

ADJACENCY_TABLE="/home/chili/chalvin/projet_long/ATLAS/ATLAS_SWORD2_filt_adjacency.tsv"

mcl $ADJACENCY_TABLE --abc -I 1.4
mcl $ADJACENCY_TABLE --abc -I 1.8
mcl $ADJACENCY_TABLE --abc -I 4.5
mcl $ADJACENCY_TABLE --abc -I 7
mcl $ADJACENCY_TABLE --abc -I 8
mcl $ADJACENCY_TABLE --abc -I 20
mcl $ADJACENCY_TABLE --abc -I 25
mcl $ADJACENCY_TABLE --abc -I 30