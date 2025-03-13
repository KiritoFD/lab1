#!/bin/bash

# Create output directory if it doesn't exist
mkdir -p /home/xy/lab1/output

# Number of files to create
read -p "How many numbered files do you want to create? " num_files

# Create the numbered files
for ((i=1; i<=$num_files; i++))
do
  echo "Creating file number $i"
  echo "This is file number $i" > "/home/xy/lab1/output/file_$i.txt"
done

echo "Created $num_files numbered files in /home/xy/lab1/output/"
