#!/usr/local/bin/python3
# -*-coding:Utf-8 -*

## Author : François-Xavier Stubbe

#
## GENETIC SCREEN
#
# Notes :
# Run the program one candidate (2 fastq files) at a time
# A "./CTRL/" folder containing a VCF (sorted, bgzip and indexed) of the CTRL strain (aligned on reference) is needed
#
# Cmd example : 
# python3 gen_screen.py -s 10E -f1 /Users/stubbf02/Fx_Stubbe/Projects/Hermand/Genetic_screen/fastq/10E_R1.fastq.gz -f2 /Users/stubbf02/Fx_Stubbe/Projects/Hermand/Genetic_screen/fastq/10E_R2.fastq.gz -g /Users/stubbf02/Fx_Stubbe/ressources/genomes/Elegans/Genomic/WBcel235.fa -bi /Users/stubbf02/Fx_Stubbe/Projects/Hermand/Genetic_screen/bwa_index/WBcel235 -q 30 -r 0.7


#################### Library import
import os
import argparse
import time
import inspect
import math

#################### Make parser

parser = argparse.ArgumentParser()

parser.add_argument("-s", "--samplename", required=True, type=str,
    help="Library name")
parser.add_argument("-f1", "--fastq1", required=True, type=str,
    help="First read file")
parser.add_argument("-f2", "--fastq2", required=True, type=str,
    help="Second read file")
parser.add_argument("-g", "--genomefasta", required=True, type=str,
    help="Full path to the genome in FASTA")
parser.add_argument("-bi", "--bwaindex", required=True, type=str,
    help="Full path to the bwa index")
parser.add_argument("-q", "--qScore", required=True, type=str,
    help="Second read file")
parser.add_argument("-r", "--ratio", required=True, type=float,
    help="Ref to Alt read ratio")

args = parser.parse_args()

#################### Catching files and variables

#Create the analysis file
os.mkdir(args.samplename)

#Files for VCF filtering
infile = "./VCF_isec/" + args.samplename + ".outann.vcf"
outfile = "./VCF_isec/" + args.samplename + "_candidates.tx"

#Catching list for ratio filtering
vcf = []   

#################### Alignment

#Trimming reads
command = "cutadapt --cores=4 -q " + args.qScore + " -o ./" + args.samplename + "/filtered_" + args.samplename + "_1" + ".fastq " + args.fastq1
os.system(command)
command = "cutadapt --cores=4 -q " + args.qScore + " -o ./" + args.samplename + "/filtered_" + args.samplename + "_2" + ".fastq " + args.fastq2
os.system(command)

#Mapping
command = "bwa mem -t 4 " + args.bwaindex + " ./" + args.samplename + "/filtered_" + args.samplename + "_1" + ".fastq" + " ./" + args.samplename +  "/filtered_" + args.samplename + "_2" + ".fastq" + " > " + "./" + args.samplename + "/aligned_" + args.samplename + ".sam"
os.system(command)

#remove unnecessary files (takes a lot of space)
os.remove(" ./" + args.samplename + "/filtered_" + args.samplename + "_1" + ".fastq")
os.remove(" ./" + args.samplename + "/filtered_" + args.samplename + "_2" + ".fastq")


#Make the BAM file
command = "samtools view -bS -T " + args.genomefasta + " ./" + args.samplename + "/aligned_" + args.samplename + ".sam > " + "./" + args.samplename + "/" + args.samplename + ".bam"
os.system(command)

#Remove the sam file
os.remove(" ./" + args.samplename + "/aligned_" + args.samplename + ".sam")

#Sort the Bam file by Chromosome position
command = "samtools sort ./" + args.samplename + "/aligned_" + args.samplename + ".bam -o ./Analysis/sorted_" + args.samplename + ".bam"
os.system(command)

#Remove the bam file
os.remove("./" + args.samplename + "/" + args.samplename + ".bam")

#Remove duplicates
command = "samtools rmdup -s ./" + samplename + "/sorted_" + args.samplename + ".bam " + "./" + samplename + "/sorted_rmdup_" + args.samplename + ".bam"
os.system(command)

#Remove sorted bam file
os.remove("./" + args.samplename + "/sorted_" + args.samplename + ".bam")

#Make summary graphs <- problem here
command = "qualimap bamqc -outdir ./" +args.samplename + "/" + args.samplename + "_Quality" +  " -bam ./" + args.samplename + "/sorted_rmdup_" + args.samplename + ".bam" + " -c -nt 4 --java-mem-size=8G"
os.system(command)

#Prepare the list for the Variant calling
with open("./" + args.samplename + "/" + args.samplename + "_list.txt" , "w") as f :
    f.write("/Users/stubbf02/Fx_Stubbe/Projects/Hermand/Genetic_screen/" + args.samplename + "/sorted_rmdup_" + args.samplename + ".bam")

#variant calling
command = "bcftools mpileup -C 50 -f " + args.genomefasta + " -O v -o ./" + args.samplename +  "/" + args.samplename + ".vcf -b ./" + args.samplename + "/" + args.samplename + "_list.txt"
os.system(command)

#################### Formatting and intersecting VCF files

#
# The idea is to get positions in the suppressors that are not in the CTRL
#

#VCF compression (bgzip required for nex step)
commad = "bgzip ./" + args.samplename + "/" + args.samplename + ".vcf"
os.system(command)

#Make vcf index
command = "bcftools index ./" + args.samplename + "/" + args.samplename + ".vcf.gz"
os.system(command)

#Get record presednt in file 1 (Supressors) but not in file 2 (CTRL)
command = "bcftools isec -C -w1 " + "./" + args.samplename + "/" + args.samplename + ".vcf.gz ./CTRL/CTRL.vcf.gz  > ./VCF_isec/" + args.samplename + ".vcf" 
os.system(command)

#Anotation of the isec VCF

command = "java -Xmx8G -jar /Users/stubbf02/Fx_Stubbe/ressources/tools/snpEff/snpEff/snpEff.jar WBcel235.75 ./VCF_isec/" + args.samplename + ".vcf > ./VCF_isec/" + args.samplename + ".outann.vcf"
os.system = command()

#################### VCF filtering

#
# The idea is to extract position in which a thresehold of mutated reads is met
# Then, annotated "MODIFIERS" are removed
#

#Filter VCF
with open(infile , "r") as b:
    with open(args.output , "w") as filehandle :
        for line in b:
            if not line.startswith("#"):
                line = line.rstrip().split("\t")
                if line[4] != "<*>" :
                    chromosome = line[0]
                    position = line[1]
                    for i in line[7].split(";"):
                        #Check the Ratio
                        if i.startswith("I16"):
                            d = i.split("=")[1].split(",")
                            refnumber = sum([int(k) for k in d[0:2]])
                            altnumber = sum([int(k) for k in d[2:4]])
                            ratio = altnumber/(altnumber+refnumber)
                    if ratio >= args.ratio:
                        refnuc = line[3]
                        altnuc = line[4].split(",")
                        #Get annotations that are not MODIFIER (so High, Medium or Low)
                        if i.startswith("ANN"):
                            annotation = i.split("=")[1].split(",")
                            for x in annotation:
                                k = x.split("|")
                                if k[2] != "MODIFIER" :
                                    altnuc = k[0]
                                    annota = k[1]
                                    purative_impact = k[2]
                                    gene_name = k[3]
                                    gene_ID = k[4]
                                    feature_type = k[5]
                                    feature_ID = k[6]
                                    transcript_biotype = k[7]
                                    rank = k[8]
                                    HGVS = k[9]
                                    HGVSp = k[10]
                                    temp = [chromosome , position , refnumber , altnumber , ratio , refnuc, altnuc, annota, purative_impact, gene_name, gene_ID, feature_type, feature_ID, transcript_biotype, rank, HGVS, HGVSp]
                                    vcf.append(temp)

#Write the output
header = ["chromosome" , "position" , "refnumber" , "altnumber" , "ratio" , "refnuc", "altnuc", "annota", "purative_impact", "gene_name", "gene_ID", "feature_type", "feature_ID", "transcript_biotype", "rank", "HGVS", "HGVSp"]
with open(outfile , "w") as filehandle :
    for columns in header: filehandle.write(columns + "\t")
    filehandle.write("\n")
    for line in vcf:
        for element in line:
            filehandle.write(str(element) + "\t")
        filehandle.write("\n")

