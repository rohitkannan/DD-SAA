#!/bin/bash

# Script for running the kNN-based reweighted SAA approach of Bertsimas and Kallus
julia-0.6.4/bin/julia main_bk.jl $1 $2 $3 $4
tar -czvf case1_bk_$1_$2_$3_$4.tar.gz case1_bk