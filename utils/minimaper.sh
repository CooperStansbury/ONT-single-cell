#!/bin/bash


FASTQ=${1}
REF=${2}
THREADS=${3}
OUTPUT=${4}

echo "FASTQ:"  "${FASTQ}"
echo "REFERENCE:" "${REF}"
echo "THREADS:" "${THREADS}"
echo "OUTPUT:" "${OUTPUT}"

minimap2 -ax splice -uf --secondary=no --MD \
-t "${THREADS}" "${REF}" "${FASTQ}" \
| samtools sort -@ "${THREADS}" -O bam -o "${OUTPUT}"

