#!/bin/bash
# organizepdfs.sh
# Organizes a folder of SPADE output PDFs into subfolders for each parameter
#
# Instructions:  Copy this script to the "pdf" folder in a SPADE output folder, and run it.
#
# Version 0.2.0
# Erin Simonds
# March 20, 2011

(
IFS=$'\n'
filenames=(`ls -1|grep .*median.*\.pdf$`)

echo "=== Unique parameters: ==="
uniqueparams=(`printf "%s\n" "${filenames[@]}" | sed 's/.*\.anno\.Rsave\.medians//g' | sed 's/\.pdf//g' |sort -u`)
printf "%s\n" "${uniqueparams[@]}"
paramstotal=(`printf "%s\n" "${uniqueparams[@]}"|wc -l`)
echo "Total unique parameters: $paramstotal"

echo "=== Unique FCS files: ==="
uniquefcs=(`printf "%s\n" "${filenames[@]}" | sed 's/\.fcs.*//g' | sort -u`)
printf "%s\n" "${uniquefcs[@]}"
fcstotal=(`printf "%s\n" "${uniquefcs[@]}"|wc -l`)
echo "Total unique FCS files: $fcstotal"

echo "=== Making folder for each parameter: ==="
mkdir ../median_pdf_organized_by_parameter
COUNT="0"
for i in ${uniqueparams[*]}; do
    mkdir -v ../median_pdf_organized_by_parameter/$i
    let COUNT++
done

echo "=== Copying each PDF to the appropriate parameter folder ==="
COUNT="0"
for i in ${uniqueparams[*]}; do
    cp -v *$i*\.pdf ../median_pdf_organized_by_parameter/$i
    let COUNT++
done

echo "=== Making folder for each FCS file: ==="
mkdir ../median_pdf_organized_by_FCS
COUNT="0"
for i in ${uniquefcs[*]}; do
    mkdir -v ../median_pdf_organized_by_FCS/$i
    let COUNT++
done

echo "=== Copying each PDF to the appropriate FCS file folder ==="
COUNT="0"
for i in ${uniquefcs[*]}; do
    cp -v $i*\.pdf ../median_pdf_organized_by_FCS/$i
    let COUNT++
done

echo "=== Finished ==="
)
