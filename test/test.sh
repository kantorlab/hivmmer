#!/bin/bash
hivmmer --id 5VM.pol --fq1 5VM_1.fastq --fq2 5VM_2.fastq --ref pol.hmm
hivmmer --id 5VM.pol --fq1 5VM_1.fastq --fq2 5VM_2.fastq --ref pol.hmm --allow-stop-codons
hivmmer --id 5VM.int --fq1 5VM_1.fastq --fq2 5VM_2.fastq --ref int.hmm --region int
