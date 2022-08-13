#!/bin/bash

# Script for running the FI-SAA instances
julia-0.6.4/bin/julia main_saa.jl $1 $2 $3 $4
tar -czvf case1_saa_$1_$2_$3_$4.tar.gz case1_saa