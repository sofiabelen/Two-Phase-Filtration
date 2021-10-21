#!/bin/bash

time julia main.jl
sxiv ../img/2phase-filtration-density.png
# sxiv ../thesis/img/img1.png
cd ../plot
python plot.py
# sxiv ../img/flows.png
