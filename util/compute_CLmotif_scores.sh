#!/bin/bash
set -e
# script to compute position-wise CL-motif scores

REF="$1"            # references sequence
BAM="$2"            # BAM file containg preprocessed aligned iCLIP/eCLIP reads
MOTIFS_XML="$3"     # output from DREME
MOTIFS_TXT="$4"
TEMP_DIR="$(mktemp -d)"       # dir to store intermediate files
BED_OUT="$5"        # BED output file containing motifs ID (column 4) and corresponding fimo score (column 5)
#BDW=$6


if [ "$#" -ne 5 ]; then
    echo "ERROR: Illegal number of arguments."
    echo "Should be 5: reference sequence, BAM file, CL-motifs.xml, CL-motifs.txt, path to temp. directory, BED outputfile."
    echo "Given arguments:"
    echo "Reference file name: " "$REF"
    echo "BAM file name: " "$BAM"
    echo "Motifs.xml file name: " "$MOTIFS_XML"
    echo "Motifs.txt file name: " "$MOTIFS_TXT"
    echo "Output file name: " "$BED_OUT"
    exit 1
fi

if [ ! -d "$TEMPDIR"]; then
  echo "Could not create temporary directory" >&2
  exit 1
else
  echo "Using temporary directory '$TEMP_DIR'"
fi

bedtools="bedtools"
if [ "$BEDTOOLS"x != "x" ]; then
  bedtools="$BEDTOOLS"
  if ! "$bedtools" --version >/dev/null 2>&1; then
    echo "Specified BEDTOOLS='$bedtools' not found."
    exit 1
  fi
else
  if ! "$bedtools" --version >/dev/null 2>&1; then
    echo "Cannot find bedtools. Specify bedtools in BEDTOOLS environment variable"
    exit 1
  fi
fi

fimo="fimo"
if [ "$FIMO"x != "x" ]; then
  fimo="$FIMO"
  if ! "$fimo" --version >/dev/null 2>&1; then
    echo "Specified FIMO='$fimo' not found."
    exit 1
  fi
else
  if ! "$fimo" --version >/dev/null 2>&1; then
    echo "Cannot find fimo. Specify fimo in FIMO environment variable"
    exit 1
  fi
fi

if [ "$WINEXTRACT"x != "x" ]; then
  if [ ! -e "$WINEXTRACT" ]; then
    echo "Specified WINEXTRACT='$WINEXTRACT' not found."
    exit 1
  fi
else
  echo "Please specify winextract in WINEXTRACT environment variable."
  exit 1
fi


#################################
# get covered regions in BAM file
# (in order to avoid running FIMO on whole genome)

# convert into bed file
"$bedtools" bamtobed -i "$BAM" > "$TEMP_DIR/alignments.bam.bed"
# only use read start positions
awk 'BEGIN{OFS=FS="\t"} $6 == "+" {print $1, ($2-1), $2, ".", ".", $6}; $6 == "-" {print $1, $3, ($3+1), ".", ".", $6} ' "$TEMP_DIR/alignments.bam.bed" > "$TEMP_DIR/alignments.bam.beginPos.bed"
# get all positions with at least one read starting
sort -k1,1 -k2,2n -k6,6 "$TEMP_DIR/alignments.bam.beginPos.bed" | uniq -c | awk 'BEGIN{OFS="\t"} { print $2,$3,$4,".", $1, $7}' > "$TEMP_DIR/alignments.bam.readStartCount.bed"

# expand each position by 3 bandwidths
w=150  #$(($BDW*3))
awk -v var=$w 'BEGIN{FS=OFS="\t"} ($2-var) >= 0 {print $1, ($2-var), ($3+var), $4, $5, $6}; ($2-var) < 0 {print $1, 0, ($3+var), $4, $5, $6};' "$TEMP_DIR/alignments.bam.readStartCount.bed" > "$TEMP_DIR/alignments.bam.expanded.bed"
sort -k1,1 -k2,2n "$TEMP_DIR/alignments.bam.expanded.bed" > "$TEMP_DIR/alignments.bam.expanded.so.bed"
# merge overlapping regions
"$bedtools" merge -s -c 6 -o distinct,distinct,distinct -i  "$TEMP_DIR/alignments.bam.expanded.so.bed" > "$TEMP_DIR/alignments.bam.regions.bed"
awk 'BEGIN{FS="\t"; OFS="\t"} {print $1, $2, $3, ".", 0, $4};' "$TEMP_DIR/alignments.bam.regions.bed" > "$TEMP_DIR/alignments.bam.regions2.bed"


######################################
# extract sequences of covered regions
"$WINEXTRACT" -g "$REF" -c "$TEMP_DIR/alignments.bam.regions2.bed" -o "$TEMP_DIR/alignment.covered_regions.sequences.fasta" -uow

######################################
# find FIMO matches
# version >= 5.0: fimo.tsv: motif_id    motif_alt_id   sequence_name   ...
# version >= 4.11.4: fimo.txt: # motif_id    motif_alt_id    sequence_name   ...
# version 4.11.3: fimo.txt: #pattern name   sequence_name    ...
echo "Run FIMO ..."
"$fimo" -oc "$TEMP_DIR/FIMO_CL_MOTIFS" --verbosity 2 --norc --motif-pseudo 0.1 --max-stored-scores 50000000 --thresh 0.01 "$MOTIFS_XML" "$TEMP_DIR/alignment.covered_regions.sequences.fasta"
if head -n1 "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.tsv" 2>/dev/null | fgrep "motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence" >/dev/null 2>&1; then
    #echo "...parse fimo output version >= 5.0"
    awk 'BEGIN{FS=OFS="\t"} $1 !~ "motif_id" && $7 > 0 {print $3, ($4-1), ($5-1), $1, $7, $6};' "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.tsv" > "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.tmp.txt"
elif head -n1 "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.txt" 2>/dev/null | fgrep "# motif_id	motif_alt_id	sequence_name	start	stop	strand	score	p-value	q-value	matched_sequence" >/dev/null 2>&1; then
    #echo "...parse fimo output version >= 4.11.4"
    awk 'BEGIN{FS=OFS="\t"} $1 !~ "# motif_id" && $7 > 0 {print $3, ($4-1), ($5-1), $1, $7, $6};' "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.txt" > "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.tmp.txt"
elif head -n1 "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.txt" 2>/dev/null | fgrep "pattern name	sequence name	start	stop	strand	score	p-value	q-value	matched sequence" >/dev/null 2>&1; then
    #echo "...parse fimo output old" 
    awk 'BEGIN{FS=OFS="\t"} $1 !~ "#pattern" && $6 > 0 {print $2, ($3-1), ($4-1), $1, $6, $5};' "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.txt" > "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.tmp.txt"
else
    echo "ERROR: Fimo output format not as expected (tested for meme 4.11.3 - 5.0.1)."
    exit 1
fi


# convert back to original positions
awk 'BEGIN{FS="\t|_"; OFS="\t"} $5 == "F" {print $2, ($3+$6), ($3+$7), $8, $9, "+"}; $5 == "R" {print $2, ($4-$7-1), ($4-$6-1), $8, $9, "-"};' "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.tmp.txt" > "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.origPos.txt"

# get single nucleotide resolution
awk 'BEGIN{FS=OFS="\t"}{chr=$1; beginPos=$2; endPos=$3; motif=$4; score=$5; strand=$6; for(i=beginPos; i <= endPos; i+=1){print chr"\t"i"\t"(i+1)"\t"motif"\t"score"\t"strand} }' "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.origPos.txt" > "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.origPos.sn.txt"

# keep only highest motif score for each position
sort -k1,1 -k2,2n -k5,5gr "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.origPos.sn.txt" > "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.origPos.sn.so.txt"
sort -k1,1 -k2,2n -s -u "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.origPos.sn.so.txt" > "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.origPos.sn.u.txt"


######################################
# replace motifs by its ID
# get DREME motifs
awk 'BEGIN{FS="\t| "; OFS="\t"} $1 == "MOTIF" {print $1, $2};' "$MOTIFS_TXT" > "$TEMP_DIR/dreme_motifs.txt"
# add id
awk -v n=$(wc -l < "$TEMP_DIR/dreme_motifs.txt") 'BEGIN{FS=OFS="\t"} {print $2, NR};' "$TEMP_DIR/dreme_motifs.txt" > "$TEMP_DIR/dreme_motifs.id.txt"
# replace motif with its id
awk 'BEGIN{FS=OFS="\t"} NR==FNR{a[$1]=$2} NR>FNR{$4=a[$4];print};' "$TEMP_DIR/dreme_motifs.id.txt" "$TEMP_DIR/FIMO_CL_MOTIFS/fimo.origPos.sn.u.txt" > "$BED_OUT"


rm -f "$TEMP_DIR/alignments.bam.bed"
rm -f "$TEMP_DIR/alignments.bam.beginPos.bed"
rm -f "$TEMP_DIR/alignments.bam.readStartCount.bed"
rm -f "$TEMP_DIR/alignments.bam.expanded.bed"
rm -f "$TEMP_DIR/alignments.bam.expanded.so.bed"
rm -f "$TEMP_DIR/alignments.bam.regions.bed"
rm -f "$TEMP_DIR/alignments.bam.regions2.bed"
rm -f "$TEMP_DIR/alignment.covered_regions.sequences.fasta"
rm -rf "$TEMP_DIR/FIMO_CL_MOTIFS"
rm -f "$TEMP_DIR/dreme_motifs.txt"
rm -f "$TEMP_DIR/dreme_motifs.id.txt"

rmdir "$TEMP_DIR"

echo "... Finished CL-motif score computation."
