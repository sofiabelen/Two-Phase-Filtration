#!/bin/bash

# rm ../dump/mass.txt
time julia main.jl
# sxiv ../img/2phase-filtration-density.png
cd ../plot
# python mass.py
# sxiv ../img/total_flow.png
