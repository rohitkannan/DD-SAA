#!/bin/bash

# Script for running the ER-SAA instances
julia-0.6.4/bin/julia main_ersaa.jl $1 $2 $3 $4
tar -czvf case1_ersaa_$1_$2_$3_$4.tar.gz case1_ersaa