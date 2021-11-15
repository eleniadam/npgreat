#!/bin/bash

#########################
#########################
# NPGREAT - ORIENTATION #
#########################
#########################
# The orientation is anchored by the telomeric nanopore reads, whose orientation is known a priori, 
# given that they always end in the 5’- (TTAGGG)n -3’ telomere tract. 
# Overlapping nanopore reads are aligned and oriented relative to the telomeric nanopore reads and to each other. 
# REXTAL sequence contigs are then oriented relative to the repeat-masked nanopore sequence reads.

###################
# Notification
echo "*** The orientation step begins... ***"

# Create folder to work inside
mkdir orientation_step
cd orientation_step

# Other folders
mkdir log_orientation
mkdir input_orientation
mkdir input_orientation/nano_tel
mkdir input_orientation/nano_subtel
mkdir input_orientation/rextal_contigs

# Input data
# Remove everything after first space in sequence ID
awk '/^>/ {$0=$1} 1' $1 > input_orientation/nano_tel/cands.fa
awk '/^>/ {$0=$1} 1' $2 > input_orientation/nano_subtel/others.fa
awk '/^>/ {$0=$1} 1' $3 > input_orientation/rextal_contigs/rextal_contigs.fa
# Extract sequences
cat input_orientation/nano_tel/cands.fa | awk '{ if (substr($0, 1, 1)==">") {filename=("input_orientation/nano_tel/" substr($0,2) ".fa")} print $0 > filename }'
cat input_orientation/nano_subtel/others.fa | awk '{ if (substr($0, 1, 1)==">") {filename=("input_orientation/nano_subtel/" substr($0,2) ".fa")} print $0 > filename }'
cat input_orientation/rextal_contigs/rextal_contigs.fa | awk '{ if (substr($0, 1, 1)==">") {filename=("input_orientation/rextal_contigs/" substr($0,2) ".fa")} print $0 > filename }'

###################
# Telomeric nanos
# Orient telomeric nanos according to p or q-end, i.e. the subtelomeric region name (4th argument)
if [[ $4 =~ "q" ]]
then
   # Keep as is
   cp input_orientation/nano_tel/cands.fa cands_or.fa
elif [[ $4 =~ "p" ]]
then
   # Reverse-Complement
   seqtk seq -r input_orientation/nano_tel/cands.fa > cands_or.fa
else
   echo "Error: The name of the subtelomeric region must specify whether it is q or p-end. For example: 10p"
fi

###################
# Subtelomeric nanos
# Move file to current folder
mv input_orientation/nano_subtel/others.fa .
# Repeat mask subtelomeric nanos
$5 -e rmblast -species human others.fa
mv others.fa.masked others_rm.fa

# Tandem repeat mask subtelomeric nanos
$6 others_rm.fa 2 7 7 80 10 50 500 -f -d -m -h
mv others_rm.fa.2.7.7.80.10.50.500.mask others_rm_trf.fa

###################
# Blastn masked subtel nanos (subject) with the oriented tel nanos (query) - pid=75
blastn -subject others_rm_trf.fa -query cands_or.fa -perc_identity 75 -outfmt "7 sstrand sseqid length" -max_hsps 1 >> output_nanos.blast
sed '/^#/ d' output_nanos.blast > strands_others_n.txt 
sort -k2,2 strands_others_n.txt > strands_others_id_sorted_n.txt

# Select orientation based on the longest alignment or the majority (2nd option)
sort -k2,2 -r -n -k3,3 strands_others_id_sorted_n.txt > strands_others_length_id_sorted_n.txt 
sort -u -k2,2 strands_others_length_id_sorted_n.txt > strands_others_length_id_sorted_n_uniq.txt

###################
# Orient subtel nanos (based on the blast output strand)
while IFS=$'\t' read -r -a line; do
    # Keep as is
    if [[ ${line[0]} == plus* ]] ; then
         cat input_orientation/nano_subtel/${line[1]}.fa >> oriented_nanos.fa
    fi
    # Reverse-Complement
    if [[ ${line[0]} == minus* ]] ; then
         seqtk seq -r input_orientation/nano_subtel/${line[1]}.fa >> oriented_nanos.fa
    fi
done < "strands_others_length_id_sorted_n_uniq.txt"

###################
# Get file of all oriented nanos
cat cands_or.fa >> oriented_nanos.fa

# Split oriented nanos in multiple files
cat oriented_nanos.fa | awk '{ if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")} print $0 > filename }'

# Save log files
mv output_nanos.blast log_orientation/
mv strands_others_length_id_sorted_n*.txt log_orientation/
# Remove not-needed files
rm others*
rm cands*
rm *.txt

###################
# Mask oriented nanos
# Repeat mask oriented nanos
$5 -e rmblast -species human oriented_nanos.fa
mv oriented_nanos.fa.masked oriented_nanos_rm.fa

# Tandem repeat mask oriented nanos
$6 oriented_nanos_rm.fa 2 7 7 80 10 50 500 -f -d -m -h
mv oriented_nanos_rm.fa.2.7.7.80.10.50.500.mask oriented_nanos_rm_trf.fa

# Remove not-needed files
rm oriented_nanos_rm.*
rm oriented_nanos.fa.*

# Notification
echo "*** The nanopore reads have been oriented. ***"

###################
###################
# Rextal contigs
# Blastn the rextal contigs (subject) with the masked oriented nanos (query) - pid=75
blastn -subject input_orientation/rextal_contigs/rextal_contigs.fa -query oriented_nanos_rm_trf.fa -perc_identity 75 -outfmt "7 sstrand sseqid length" -max_hsps 1 > output_contigs.blast
sed '/^#/ d' output_contigs.blast > strands_others_c.txt
sort -k2,2 strands_others_c.txt > strands_others_id_sorted_c.txt

# Select orientation based on the longest alignment
sort -k2,2 -r -n -k3,3 strands_others_id_sorted_c.txt > strands_others_length_id_sorted_c.txt
sort -u -k2,2 strands_others_length_id_sorted_c.txt > strands_others_length_id_sorted_c_uniq.txt

###################
# Orient contigs (based on the blast output strand)
while IFS=$'\t' read -r -a line; do
    # Keep as is
    if [[ ${line[0]} == plus* ]] ; then
         cat input_orientation/rextal_contigs/${line[1]}.fa >> oriented_rextal_contigs.fa
    fi
    # Reverse-Complement
    if [[ ${line[0]} == minus* ]] ; then
         seqtk seq -r input_orientation/rextal_contigs/${line[1]}.fa >> oriented_rextal_contigs.fa
    fi
done < "strands_others_length_id_sorted_c_uniq.txt"

###################
# Split oriented nanos in multiple files
cat oriented_rextal_contigs.fa | awk '{ if (substr($0, 1, 1)==">") {filename=(substr($0,2) ".fa")} print $0 > filename }'

# Save log files
mv output_contigs.blast log_orientation/
mv strands_others_length_id_sorted_c*.txt log_orientation/
# Remove not-needed files
rm *.txt
rm -r input_orientation

# Notification
echo "*** The REXTAL contigs have been oriented. ***"
echo "*** The orientation step finished. ***"
echo "**********************************************"

