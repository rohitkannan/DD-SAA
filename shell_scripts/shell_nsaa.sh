#!/bin/bash

# Script for running the naive SAA instances
julia-0.6.4/bin/julia main_nsaa.jl $1 $2 $3 $4
tar -czvf case1_nsaa_$1_$2_$3_$4.tar.gz case1_nsaa