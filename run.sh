#!/bin/sh

#SBATCH -J solo
#SBATCH  -p fat
#SBATCH --nodes=1
#SBATCH --ntasks-per-node=1

python solo_gRNA.py offtarget --genome /public/home/wangycgroup/wuj/MyDb/Genome/Hsa/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --bowtie_index /public/home/wangycgroup/wuj/MyDb/Genome/Hsa/GRCh38/bowtie_idx/Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf /public/home/wangycgroup/public/00_Genome_ref/Homo_sapiens/Homo_sapiens.GRCh38.105.gtf --grnas bm.grna.fa --region rRNA_regions.bed
#python solo_gRNA.py solo --manual_ann rRNA_regions.bed --genome /public/home/wangycgroup/wuj/MyDb/Genome/Hsa/GRCh38/Homo_sapiens.GRCh38.dna.primary_assembly.fa --bowtie_index /public/home/wangycgroup/wuj/MyDb/Genome/Hsa/GRCh38/bowtie_idx/Homo_sapiens.GRCh38.dna.primary_assembly.fa --gtf /public/home/wangycgroup/public/00_Genome_ref/Homo_sapiens/Homo_sapiens.GRCh38.105.gtf --grnas input.fa
