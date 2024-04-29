#!/bin/bash

# Retreave the second best hit results and save as a text file
awk '$1==$2 {getline ; if($1!=$2) print $0}' SelfBlast > Next_Best_Hit.txt

# Get the sequence length of each alignment
awk '{print $4}' Next_Best_Hit.txt > align_length.txt

# Get the numbers of mismatches for each second best alignment
awk '{print $5}' Next_Best_Hit.txt > mismatches.txt

# Get the numbers of gap openings for each second best alignment
awk '{print $6}' Next_Best_Hit.txt > gap_openings.txt





