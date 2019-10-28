# !/bin/bash

### Maria Mikedis mikedis@wi.mit.edu

### for each analysis, update variable DIR, SAMPLE1, SAMPLE2, SAMPLE3, TPMS, EXPRESSED, PYTHDIR, FA, ANNO, BAMDIR
### This script will create directory <${ASPEAK_output_directory}/foundpeaks/merged.peaks.190622> and put files into region-specific directories


DIR="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/20180730_ASPeak/ASPeak_IgG_RNAseq_gapno0_CL.site.only_w.intergenic"



TPMS="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/20180730_3S_RNAseq_kallisto/190622_kallisto/kallisto/MM257_3S.undiff.gonia_kallisto_mm10_GRCm38_mRNA_retrogenes_no.rRNA_TPMs_190622_min.TPM.1_gene.id.only"
### file contains a list of coding transcripts, retrogenes, and ncRNAs (not rRNA) that are expressed at least at 1 TPM

EXPDIR="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/20180730_3S_RNAseq_kallisto/190622_kallisto/kallisto"
EXPFILE="mm10_GRCm38_MM257.tpm1_coding.utr3.exons"
	### Directory and base name of bed file (i.e., actual file is ${EXPFILE}.bed)
	###bed file contains 3' UTRs from genes expressed at a minimum TPM of 1; 
	### if there are multiple isoforms, use the 3' UTR of the longest isoform with a minimum TPM of 1
	### if no isoform meets cut off of 1 TPM, then use longest isoform with any number of read counts
	### to be used as background for motif analysis of 3'UTR peaks 
EXPTRANS="MM257_3S.undiff.gonia_kallisto_mm10_GRCm38_mRNA_retrogenes_no.rRNA_TPMs_190622_genes.transcripts_min.TPM.1"
	### same as above; one transcript per gene
	### contains transcript ID, gene ID, transcript TPMs for samples and mean, adn gene TPM for samples and mean

PYTHDIR="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/python_scripts_post_ASPeak_analysis"
	### directory containing the following python scripts for post-ASPeak analysis (these python scripts are called within this bash script)

RDIR="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/R_scripts_post_ASPeak_analysis"

FA="/nfs/genomes/mouse_mm10_dec_11_no_random/fasta_whole_genome/mm10.fa"
	### fasta of chromosome sequences (each chromosome is a fasta entry)

ANNO="/lab/solexa_page/maria/genome_files/kallisto/for_ASPeak_analysis/mm10_GRCm38_refGene_mRNA_retrogenes_names.txt"
	### annotation file of transcript IDs (field 1) and gene IDs (field 2)

BAMDIR="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/20180416_iCLIP_STAR"
	### directory containing bam files of aligned CLIP data

SAMPLE1="MM218_1"
SAMPLE2="MM222_1_1"
SAMPLE3="MM222_2_1"



cd ${DIR}/foundpeaks
mkdir merged.peaks.190622


### convert peak files to bed files; merge peaks so that one peak can represent multiple transcripts; link to merged.peak directory
### if no peaks meet the designated cut offs, bedtools merge will give the error message: "ERROR: Requested column 4, but database file stdin only has fields 1 - 0."
### bedtools output the strand in the wrong column so modifying bedtools output to follow bed format
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only 
	do
	mkdir "${DIR}/foundpeaks/merged.peaks.190622/${REGION}/"
	DIR2="${DIR}/foundpeaks/demultiplexed_${SAMPLE1}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${DIR2}
	awk '{FS="<"; if (NR==1) print $0; else print $1$2}' mm10_GRCm38_${REGION}.reg.peaks| awk '{FS="\t"; OFS="\t"; if (NR==1) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12,$13; if ($13 < 0.05) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12,$13}' | sort -k1,1 -k2,2n| bedtools merge -header -s -c 4,5,7,8,9,10,11,12,13 -o distinct,max,max,max,min,max,max,max,max -i stdin | awk '{FS="\t"; OFS="\t"; if (NR==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13; else print $1,$2,$3,$5,$6,$4,$7,$8,$9,$10,$11,$12,$13}' | uniq >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed;
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/
	ln -s ${DIR2}/${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed .
	DIR2="${DIR}/foundpeaks/demultiplexed_${SAMPLE2}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${DIR2}
	awk '{FS="<"; if (NR==1) print $0; else print $1$2}' mm10_GRCm38_${REGION}.reg.peaks| awk '{FS="\t"; OFS="\t"; if (NR==1) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12,$13; if ($13 < 0.05) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12,$13}' | sort -k1,1 -k2,2n| bedtools merge -header -s -c 4,5,7,8,9,10,11,12,13 -o distinct,max,max,max,min,max,max,max,max -i stdin | awk '{FS="\t"; OFS="\t"; if (NR==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13; else print $1,$2,$3,$5,$6,$4,$7,$8,$9,$10,$11,$12,$13}' | uniq >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed;
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/
	ln -s ${DIR2}/${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed .
	DIR2="${DIR}/foundpeaks/demultiplexed_${SAMPLE3}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${DIR2}
	awk '{FS="<"; if (NR==1) print $0; else print $1$2}' mm10_GRCm38_${REGION}.reg.peaks| awk '{FS="\t"; OFS="\t"; if (NR==1) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12,$13; if ($13 < 0.05) print $1,$8,$9,$3,$6,$2,$4,$5,$7,$10,$11,$12,$13}' | sort -k1,1 -k2,2n| bedtools merge -header -s -c 4,5,7,8,9,10,11,12,13 -o distinct,max,max,max,min,max,max,max,max -i stdin | awk '{FS="\t"; OFS="\t"; if (NR==1) print $1,$2,$3,$4,$5,$6,$7,$8,$9,$10,$11,$12,$13; else print $1,$2,$3,$5,$6,$4,$7,$8,$9,$10,$11,$12,$13}' | uniq >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed;
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/
	ln -s ${DIR2}/${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed .
done

### for peaks that are found in multiple "regions", remove peaks based on the following hierarchy: 3' UTR coding only > 5' UTR coding only > CDS > ncRNA > retrogenes > intron 
### for example, if a peak shows up in both 3' UTR and ncRNA, the hierarchical priortization removes it from ncRNA and leaves it in 3' UTR
### if there is only partial overlap in the peaks, the full peak will be removed

for REGION in intron_coding.only; do 
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}

	bedtools subtract -header -a ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE1}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE1}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${SAMPLE1}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| temp4;
	bedtools subtract -header -a temp4 -b ../retrogenes/${SAMPLE1}_mm10_GRCm38_retrogenes.reg.peaks.bed -s -A >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE2}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE2}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${SAMPLE2}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| temp4;
	bedtools subtract -header -a temp4 -b ../retrogenes/${SAMPLE2}_mm10_GRCm38_retrogenes.reg.peaks.bed -s -A >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE3}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE3}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${SAMPLE3}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| temp4;
	bedtools subtract -header -a temp4 -b ../retrogenes/${SAMPLE3}_mm10_GRCm38_retrogenes.reg.peaks.bed -s -A >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done

for REGION in retrogenes; do 
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}

	bedtools subtract -header -a ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE1}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE1}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${SAMPLE1}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE2}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE2}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${SAMPLE2}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE3}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE3}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| temp3;
	bedtools subtract -header -a temp3 -b ../ncRNA.no.rRNA/${SAMPLE3}_mm10_GRCm38_ncRNA.no.rRNA.reg.peaks.bed -s -A >|  ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done

for REGION in ncRNA.no.rRNA; do 
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}

	bedtools subtract -header -a ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE1}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE1}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE2}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE2}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE3}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| temp2;
	bedtools subtract -header -a temp2 -b ../cds/${SAMPLE3}_mm10_GRCm38_cds.reg.peaks.bed -s -A >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done

for REGION in cds; do 
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}

	bedtools subtract -header -a ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE1}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE2}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >|  ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

	bedtools subtract -header -a ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| temp1;
	bedtools subtract -header -a temp1 -b ../utr5_coding.only/${SAMPLE3}_mm10_GRCm38_utr5_coding.only.reg.peaks.bed -s -A >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
	rm -f temp*;

done



for REGION in utr5_coding.only; do 
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}

	bedtools subtract -header -a ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE1}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

	bedtools subtract -header -a ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE2}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

	bedtools subtract -header -a ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ../utr3_coding.only/${SAMPLE3}_mm10_GRCm38_utr3_coding.only.reg.peaks.bed -s -A >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

done

	### saving utr3_coding.only peak file with hierachry.applied in name to keep file names consistent; no utr3 peaks are removed
for REGION in utr3_coding.only; do 
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}

	cat ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

	cat ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

	cat ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;

done



### identify peaks that are present in at least 2 replicates using bedtools; one loop for regions where hierachy was applied; second loop for other regions
for REGION in cds intron_coding.only ncRNA.no.rRNA retrogenes utr3_coding.only utr5_coding.only; do
cd  ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/
CLIP1=${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
CLIP2=${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
CLIP3=${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed;
bedtools intersect -header -s -a ${CLIP2} -b ${CLIP1}| sort -k1,1 -k2,2n >| ${SAMPLE1}_${SAMPLE2}_intersect_${REGION}.bed;
bedtools intersect -header -s -a ${CLIP1} -b ${CLIP3}| sort -k1,1 -k2,2n >| ${SAMPLE1}_${SAMPLE3}_intersect_${REGION}.bed;
bedtools intersect -header -s -a ${CLIP3} -b ${CLIP2}| sort -k1,1 -k2,2n >| ${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed;
bedtools intersect -header -s -a ${SAMPLE1}_${SAMPLE2}_intersect_${REGION}.bed -b ${CLIP3}| sort -k1,1 -k2,2n >| ${SAMPLE1}_${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed;
done

for REGION in intergenic rRNA.only tRNAs; do
cd  ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/
CLIP1=${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed;
CLIP2=${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed;
CLIP3=${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed;
bedtools intersect -header -s -a ${CLIP2} -b ${CLIP1}| sort -k1,1 -k2,2n >| ${SAMPLE1}_${SAMPLE2}_intersect_${REGION}.bed;
bedtools intersect -header -s -a ${CLIP1} -b ${CLIP3}| sort -k1,1 -k2,2n >| ${SAMPLE1}_${SAMPLE3}_intersect_${REGION}.bed;
bedtools intersect -header -s -a ${CLIP3} -b ${CLIP2}| sort -k1,1 -k2,2n >| ${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed;
bedtools intersect -header -s -a ${SAMPLE1}_${SAMPLE2}_intersect_${REGION}.bed -b ${CLIP3}| sort -k1,1 -k2,2n >| ${SAMPLE1}_${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed;
done


### create a master bed file with genomic position information only for peaks that are present in at least 2 of the 3 biological replicates
### I used the following strategy instead of bedtools because bedtools merge was doing funny things when I tried to merge intersected files
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
cd  ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/
cut -f1-3,6 ${SAMPLE1}_${SAMPLE2}_intersect_${REGION}.bed | sort -k1,1 -k2,2n| uniq >| tmp1.tmp
cut -f1-3,6 ${SAMPLE1}_${SAMPLE3}_intersect_${REGION}.bed| sort -k1,1 -k2,2n | uniq >| tmp2.tmp
cut -f1-3,6 ${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed | sort -k1,1 -k2,2n| uniq  >| tmp3.tmp
cat tmp1.tmp tmp2.tmp >| tmp.tmp
cat tmp.tmp tmp3.tmp | sort -k1,1 -k2,2n | uniq | awk '{FS="\t"; print $1"\t"$2"\t"$3"\t""name""\t""score""\t"$4}'>|  ${REGION}.peaks_present.in.at.least.2.replicates_peak.info.only.bed
rm -f tmp*
done



####pull out replicated peaks from each replicates peak.bed file 
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
	REPS=${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${REGION}.peaks_present.in.at.least.2.replicates_peak.info.only.bed;
	DIR2="${DIR}/foundpeaks/demultiplexed_${SAMPLE1}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${DIR2};
	REPS=${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${REGION}.peaks_present.in.at.least.2.replicates_peak.info.only.bed;
	bedtools intersect -header -s -a ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ${REPS} >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed;
	DIR2="${DIR}/foundpeaks/demultiplexed_${SAMPLE2}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${DIR2};
	bedtools intersect -header -s -a ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ${REPS} >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed;
	DIR2="${DIR}/foundpeaks/demultiplexed_${SAMPLE3}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks";
	cd ${DIR2};
	bedtools intersect -header -s -a ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks.bed -b ${REPS} >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed;
done

### combine bare-bones master bed file with peak information to create a master bed file with more complete info

for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}	
	FILE1="${DIR}/foundpeaks/demultiplexed_${SAMPLE1}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks/${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed";
	FILE2="${DIR}/foundpeaks/demultiplexed_${SAMPLE2}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks/${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed";
	FILE3="${DIR}/foundpeaks/demultiplexed_${SAMPLE3}_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted/peaks/${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_replicated.bed";
	cat ${FILE1} ${FILE2} >| tmp1.tmp
	cat tmp1.tmp ${FILE3} | sort -k1,1 -k2,2n | bedtools merge -header -s -c 4,5,6,7,8,9,10,11,12,13 -o distinct,max,distinct,max,max,min,max,max,max,max -i stdin | uniq |awk '{FS="\t"; OFS="\t"; print $1,$2,$3,$5,$6,$7,$8,$9,$10,$11,$12,$13}' >| ${REGION}.peaks_present.in.at.least.2.replicates.bed; # uniq command gets rid of repeated headers; bedtools behaving weird and output 2 strand columns
	rm -f tmp*;
done



### create a list of genes with peaks, without header

for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only tRNAs utr3_coding.only utr5_coding.only; do
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}
	cut -f4 ${REGION}.peaks_present.in.at.least.2.replicates.bed | awk '{FS="_"; if (NR!=1) print $1}' | sort | uniq >| ${REGION}.peaks_present.in.at.least.2.replicates_gene.ID.only;
done


for REGION in retrogenes; do 
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}
	cut -f4 ${REGION}.peaks_present.in.at.least.2.replicates.bed | awk '{FS="_exon"; if (NR!=1) print $1}' | sort | uniq >| ${REGION}.peaks_present.in.at.least.2.replicates_gene.ID.only;
done


### calculate overlap in peaks among replicates; produce venn diagrams
for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}
	area1="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed | cut -f1-3,6| uniq | wc | awk '{print $1}')"
	area2="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed | cut -f1-3,6| uniq| wc | awk '{print $1}')"
	area3="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed | cut -f1-3,6| uniq| wc | awk '{print $1}')"
	n12="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${SAMPLE1}_${SAMPLE2}_intersect_${REGION}.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	n23="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	n13="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${SAMPLE1}_${SAMPLE3}_intersect_${REGION}.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	n123="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${SAMPLE1}_${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	echo "Replicate No_peaks" >| tmp.tmp
	echo "${SAMPLE1} ${area1}" >> tmp.tmp
	echo "${SAMPLE2} ${area2}">> tmp.tmp
	echo "${SAMPLE3} ${area3}" >> tmp.tmp
	echo "${SAMPLE1}_${SAMPLE2} ${n12}" >> tmp.tmp
	echo "${SAMPLE2}_${SAMPLE3} ${n23}" >> tmp.tmp
	echo "${SAMPLE1}_${SAMPLE3} ${n13}" >> tmp.tmp
	echo "${SAMPLE1}_${SAMPLE2}_${SAMPLE3} ${n123}" >> tmp.tmp
	awk '{if ($2=="") print $1"\t""0"; else print $1"\t"$2}' tmp.tmp >| ${REGION}_overlapping_peaks
	rm -f *.tmp
	Rscript ${RDIR}/make_venn_diagram.R ${REGION}
done


### retrieve sequences at replicated peaks; make a graph of nucleotides represented at peaks; graphs produced assume that the peaks are at crosslinked sites, which may not be true

for REGION in cds intron_coding.only intergenic ncRNA.no.rRNA rRNA.only retrogenes tRNAs utr3_coding.only utr5_coding.only; do
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}
	mkdir sequence_analysis
	cd sequence_analysis
	bedtools getfasta -s -fi /lab/solexa_page/maria/genome_files/mm10_chrMT.fa -bed ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/${REGION}.peaks_present.in.at.least.2.replicates.bed |grep -v ">"| sed 's/\(.\)/\1\n/g'| sed '/^$/d' | sed 's/.*/\U&/' >| ${REGION}.peaks_present.in.at.least.2.replicates_CL.nts.only ### format is one nucleotide per line, representing a crosslinked site
	Rscript /lab/solexa_page/maria/scripts/CLIP_analysis/nucleotides_at_crosslinked_sites.R ${REGION}
done




### for regions that were quantified via RNA-seq (coding transcripts, ncRNA, retrogene), limit called peaks to transcripts with at least a TPM of 1
### grep is very slow but I couldn't figure out a faster version of grep that actually worked
for REGION in cds retrogenes utr3_coding.only utr5_coding.only; do
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}
	mkdir min.TPM.1
	cd min.TPM.1
	awk 'BEGIN {FS="\t"} {print $0"_"} ' ${TPMS} >| TPMs_temp ## to make sure grep pulls out exact matches
	while read name; do
	  grep "\<$name" ../${REGION}.peaks_present.in.at.least.2.replicates.bed 
	done <TPMs_temp >|  ${REGION}.peaks_present.in.at.least.2.replicates_min.TPM.1.bed 
	rm -f TPMs_temp
done

	### ncRNA gene names are labeled with "(ncRNA)" to differentiate from coding transcripts with the same gene name
	### pull out gene names with ncRNA label and then see if only those are expressed
for REGION in ncRNA.no.rRNA; do
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}
	mkdir min.TPM.1
	cd min.TPM.1
	awk 'BEGIN {FS=" "} {if ($2=="(ncRNA)") print $1"_"} ' ${TPMS} >| TPMs_temp ## to make sure file pulls out exact match
	while read name; do
	  grep "\<$name" ../${REGION}.peaks_present.in.at.least.2.replicates.bed
	done <TPMs_temp  >| ${REGION}.peaks_present.in.at.least.2.replicates_min.TPM.1.bed 
	rm -f TPMs_temp
done




#### repeat previous analysis using only transcripts whose genes are expressed at least at TPM of 1
### create a list of genes with peaks, without header
for REGION in cds ncRNA.no.rRNA utr3_coding.only utr5_coding.only; do
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/min.TPM.1
	cut -f4 ${REGION}.peaks_present.in.at.least.2.replicates_min.TPM.1.bed | awk '{FS="_"; if (NR!=1) print $1}' | sort | uniq >| ${REGION}.peaks_present.in.at.least.2.replicates_min.TPM.1_gene.ID.only;
done


for REGION in retrogenes; do 
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/min.TPM.1
	cut -f4 ${REGION}.peaks_present.in.at.least.2.replicates_min.TPM.1.bed | awk '{FS="_exon"; if (NR!=1) print $1}' | sort | uniq >| ${REGION}.peaks_present.in.at.least.2.replicates_min.TPM.1_gene.ID.only;
done

### identify overlap peaks per replicate
### calculate overlap in peaks among replicates; produce venn diagrams
### grep is very slow but I couldn't figure out a faster version of grep that actually worked
### The total number of replicated peaks in the venn diagram is slighly greater than the number in the merged bed file because consecutive nts may be present in different 2 out of 3 datasets; these consecutive peaks are merged into 1 peak within the bed file
for REGION in cds ncRNA.no.rRNA retrogenes utr3_coding.only utr5_coding.only; do
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/min.TPM.1
	awk 'BEGIN {FS="\t"} {print $0"_"} ' ${TPMS} >| TPMs_temp
	while read name; do
	  grep "\<$name" ../${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed 
	done <TPMs_temp >| ${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied_min.TPM.1.bed
	while read name; do
	  grep "\<$name" ../${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed 
	done <TPMs_temp >| ${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied_min.TPM.1.bed
	while read name; do
	  grep "\<$name" ../${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied.bed 
	done <TPMs_temp >| ${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied_min.TPM.1.bed
	while read name; do
	  grep "\<$name" ../${SAMPLE1}_${SAMPLE2}_intersect_${REGION}.bed 
	done <TPMs_temp >| ${SAMPLE1}_${SAMPLE2}_intersect_${REGION}_min.TPM.1.bed 
	while read name; do
	  grep "\<$name" ../${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed 
	done <TPMs_temp >| ${SAMPLE2}_${SAMPLE3}_intersect_${REGION}_min.TPM.1.bed
	while read name; do
	  grep "\<$name" ../${SAMPLE1}_${SAMPLE3}_intersect_${REGION}.bed 
	done <TPMs_temp >| ${SAMPLE1}_${SAMPLE3}_intersect_${REGION}_min.TPM.1.bed
	while read name; do
	  grep "\<$name" ../${SAMPLE1}_${SAMPLE2}_${SAMPLE3}_intersect_${REGION}.bed 
	done <TPMs_temp >| ${SAMPLE1}_${SAMPLE2}_${SAMPLE3}_intersect_${REGION}_min.TPM.1.bed
	area1="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/min.TPM.1/${SAMPLE1}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied_min.TPM.1.bed | cut -f1-3,6| uniq | wc | awk '{print $1}')"
	area2="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/min.TPM.1/${SAMPLE2}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied_min.TPM.1.bed | cut -f1-3,6| uniq| wc | awk '{print $1}')"
	area3="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/min.TPM.1/${SAMPLE3}_mm10_GRCm38_${REGION}.reg.peaks_hierarchy.applied_min.TPM.1.bed | cut -f1-3,6| uniq| wc | awk '{print $1}')"
	n12="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/min.TPM.1/${SAMPLE1}_${SAMPLE2}_intersect_${REGION}_min.TPM.1.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	n23="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/min.TPM.1/${SAMPLE2}_${SAMPLE3}_intersect_${REGION}_min.TPM.1.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	n13="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/min.TPM.1/${SAMPLE1}_${SAMPLE3}_intersect_${REGION}_min.TPM.1.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	n123="$(grep -v "#" ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/min.TPM.1/${SAMPLE1}_${SAMPLE2}_${SAMPLE3}_intersect_${REGION}_min.TPM.1.bed | cut -f1-3,6| sort -k1,1 -k2,2n | uniq| wc | awk '{print $1}')"
	echo "Replicate No_peaks" >| tmp.tmp
	echo "${SAMPLE1} ${area1}" >> tmp.tmp
	echo "${SAMPLE2} ${area2}">> tmp.tmp
	echo "${SAMPLE3} ${area3}" >> tmp.tmp
	echo "${SAMPLE1}_${SAMPLE2} ${n12}" >> tmp.tmp
	echo "${SAMPLE2}_${SAMPLE3} ${n23}" >> tmp.tmp
	echo "${SAMPLE1}_${SAMPLE3} ${n13}" >> tmp.tmp
	echo "${SAMPLE1}_${SAMPLE2}_${SAMPLE3} ${n123}" >> tmp.tmp
	awk '{if ($2=="") print $1"\t""0"; else print $1"\t"$2}' tmp.tmp >| ${REGION}_overlapping_peaks
	rm -f *.tmp
	Rscript ${RDIR}/make_venn_diagram.R ${REGION}
done






### retrieve sequences at replicated peaks; make a graph of nucleotides represented at peaks

for REGION in cds ncRNA.no.rRNA retrogenes utr3_coding.only utr5_coding.only; do
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/min.TPM.1
	mkdir sequence_analysis
	cd sequence_analysis
	bedtools getfasta -s -fi /lab/solexa_page/maria/genome_files/mm10_chrMT.fa -bed ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/min.TPM.1/${REGION}.peaks_present.in.at.least.2.replicates_min.TPM.1.bed |grep -v ">"| sed 's/\(.\)/\1\n/g'| sed '/^$/d' | sed 's/.*/\U&/' >| ${REGION}.peaks_present.in.at.least.2.replicates_CL.nts.only ### format is one nucleotide per line, representing a crosslinked site
	Rscript /lab/solexa_page/maria/scripts/CLIP_analysis/nucleotides_at_crosslinked_sites.R ${REGION}
done




### motif analysis of replicated peaks in 3' UTR only, using homer
### background: 3' UTRs from genes expressed at a minimum TPM of 1
### if there are multiple isoforms, use the 3' UTR of the most abundantly expressed isoform
	#### step 1: create background bed file for HOMER analysis
	### get fasta from bed file
	### get a file of listed RefSeq IDs with length of 3' UTR
cd ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis
mkdir background_transcriptome
cd background_transcriptome
ln -s ${EXPDIR}/${EXPFILE}.bed . 

	### bsub -q 18 sends job to ubuntu cluster 18 nodes
bsub -q 18 -K Rscript ${RDIR}/extract_3UTR_sequences.R 3UTRs sleep 10 &
wait

### remove line breaks in fasta file that was just downloaded
 awk '!/^>/ { printf "%s", $0; n = "\n" } 
/^>/ { print n $0; n = "" }
END { printf "%s", n }
' 3UTRs.fa >| temp.fa
rm -f 3UTRs.fa
mv temp.fa 3UTRs.fa


ln -s ${EXPDIR}/MM257_3S.undiff.gonia_kallisto_mm10_GRCm38_mRNA_retrogenes_no.rRNA_TPMs_190622_genes.transcripts_min.TPM.1_most.expressed.coding.isoform.per.gene . 
TRANSCRIPTS="MM257_3S.undiff.gonia_kallisto_mm10_GRCm38_mRNA_retrogenes_no.rRNA_TPMs_190622_genes.transcripts_min.TPM.1_most.expressed.coding.isoform.per.gene"



### Create fasta with 3'UTR for most expressed isoform per gene
awk 'BEGIN {FS="\t"} {print ">"$2"+"}' ${TRANSCRIPTS} | sort | uniq >| temp.Refseq.ids 
awk 'BEGIN {FS="\t"} {print $1"+"}' 3UTRs.fa | grep -A 1 -Ff temp.Refseq.ids | awk 'BEGIN {FS="+"} {print $1}' | sed '/--/d' >| 3UTRs_most.expressed.coding.isoform.per.gene.fa    ### sed removes lines that are marked as "--"



### get bed file of 3' UTRs from the most expressed isoforms per gene
awk 'BEGIN {FS="\t"} {print $2"_"}' ${TRANSCRIPTS} | sort | uniq >| temp.Refseq.ids2
grep -Ff temp.Refseq.ids2 ${EXPFILE}.bed >| ${EXPFILE}_most.expressed.coding.isoform.per.gene.bed 


### get bed file of 3' UTRs from the most expressed isoforms per gene that are DAZL targets only 
awk 'BEGIN {FS="\t"} {print $1"_"}' utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_gene.ID.only >| temp.gene.id
awk 'BEGIN {FS="\t"} {print $1"_\t"$0}' ${TRANSCRIPTS}	| grep -Ff temp.gene.id | awk 'BEGIN {FS="\t"} {print $3"_"}' >|  temp.Refseq.ids3
grep -Ff temp.Refseq.ids3 ${EXPFILE}.bed >| ${EXPFILE}_most.expressed.coding.isoform.per.gene_DAZL.bound.bed 

###  get fa file of 3' UTRs from the most expressed isoforms per gene that are DAZL targets only 
awk 'BEGIN {FS="\t"} {print $1"_\t"$0}' ${TRANSCRIPTS}	| grep -Ff temp.gene.id | awk 'BEGIN {FS="\t"} {print ">"$3"+"}' >|  temp.Refseq.ids4
awk 'BEGIN {FS="\t"} {print $1"+"}' 3UTRs.fa | grep -A 1 -Ff temp.Refseq.ids4 | awk 'BEGIN {FS="+"} {print $1}' | sed '/--/d' >| 3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound.fa    ### sed removes lines that are marked as "--"




rm -f temp*

cd ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis
mkdir homer_output_background.exp.transcripts


	### and finally, motif analysis using HOMER; homer must be run from its own directory and cannot be run on cluster
cd /lab/solexa_page/maria/homer/bin/

INPUT="${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1.bed"
GENOME="mm10"
OUTPUT="${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/homer_output_background.exp.transcripts"
BKGD="${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/background_transcriptome/${EXPFILE}_most.expressed.coding.isoform.per.gene.bed"
/lab/solexa_page/maria/homer/bin/findMotifsGenome.pl ${INPUT} ${GENOME} ${OUTPUT} -bg ${BKGD} -size 21 -S 10 -len 3,4,5,6, -rna -preparse -preparsedDir ${OUTPUT}










#### motif analysis via MEME

cd ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis
mkdir meme_analysis
cd meme_analysis

ln -s ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1.bed . 

### expand bed file by 20 nucleotides on each side of crosslinked sites; extract as fasta sequences
	### "-q 18" send job to ubuntu cluster 18 nodes
### R script requries a list of transcripts from which I want to extract sequences
cut -f4 ../background_transcriptome/mm10_GRCm38_MM257.tpm1_coding.utr3.exons_most.expressed.coding.isoform.per.gene.bed | awk 'BEGIN {FS="_"} {print $1"_"$2}' | sort | uniq >| temp.Refseq.ids
bsub -q 18 -K Rscript ${RDIR}/obtain_3UTRseq_from_bed_190220.R utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1.bed temp.Refseq.ids 20 temp sleep 10 &
wait

	### shorten header names so that meme can deal with them more easily; convert to RNA sequence
awk 'BEGIN {FS="\n"} {if (substr($1,1,1) == ">") print ">"NR; else print $0;} ' temp.fa | sed 's/T/U/g' >| utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_expanded20.rna_header.modified.fa

### use background fasta file that was already create in background directory; convert to RNA sequence
sed 's/T/U/g' ../background_transcriptome/3UTRs_most.expressed.coding.isoform.per.gene.fa >| 3UTRs_most.expressed.coding.isoform.per.gene.rna.fa



BKGDFASTA="3UTRs_most.expressed.coding.isoform.per.gene.rna.fa"
bsub -q 18 -K /usr/local/meme/bin/fasta-get-markov -rna -m 0  ${BKGDFASTA} ${EXPFILE}_most.expressed.coding.isoform.per.gene_markov.model.0.order.txt sleep 10 &
bsub -q 18 -K /usr/local/meme/bin/fasta-get-markov -rna -m 1  ${BKGDFASTA} ${EXPFILE}_most.expressed.coding.isoform.per.gene_markov.model.1.order.txt sleep 10 &
wait

### motif analysis; note that consecutive crosslinked sites are merged so will only be analyzed once
### background: all expressed 3' UTRs (gene TPM>=1); if more than one isoform is expressed, most robustly expressed isoform used
DAZLBOUND="utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_expanded20.rna_header.modified.fa"
bsub -n 16 "/usr/local/meme/bin/meme ${DAZLBOUND} -rna -oc output_bkgd.all.expressed.markov.0.order_maxw6 -mod oops -nmotifs 6 -minw 3 -maxw 6 -maxsize 2000000 -p 16 -bfile ${EXPFILE}_most.expressed.coding.isoform.per.gene_markov.model.0.order.txt"
bsub -n 16 "/usr/local/meme/bin/meme ${DAZLBOUND} -rna -oc output_bkgd.all.expressed.markov.1.order_maxw6 -mod oops -nmotifs 6 -minw 3 -maxw 6 -maxsize 2000000 -p 16 -bfile ${EXPFILE}_most.expressed.coding.isoform.per.gene_markov.model.1.order.txt"

rm -f temp*

###################
### metaplot analysis around crosslinked site
###################

cd ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis
mkdir metaplot_around_CL_nt
cd metaplot_around_CL_nt

### convert bam files containing crosslinked position to bed files
for FILE in demultiplexed_MM218_1_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted demultiplexed_MM222_1_1_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted demultiplexed_MM222_2_1_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted; do 
bedtools bamtobed -i ${BAMDIR}/${FILE}.bam >| ${FILE}.bed 
done 

### combine crosslinked sites from 3 biological replicates
for FILE in demultiplexed_MM218_1_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted demultiplexed_MM222_1_1_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted demultiplexed_MM222_2_1_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted; do 
awk 'BEGIN {FS="\t"} {print $1"\t"$2"\t"$3"\t""n/a""\t"$5"\t"$6}' ${FILE}.bed >> CL.site.only_from_bam_files.bed
sort -k1,1 -k2,2g CL.site.only_from_bam_files.bed | uniq >| tmp 
mv -f tmp CL.site.only_from_bam_files.bed
done

rm -f *DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted.bed

### pull out replicated crosslink sites, one nt per line
DIR2="${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/"
bedtools intersect -a ${DIR2}/utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1.bed -b CL.site.only_from_bam_files.bed >| utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line.bed


### expand bed file by 17 nucleotides on each side of crosslinked sites; extract as fasta sequences
	### "-q 18" send job to ubuntu cluster 18 nodes
### R script requries a list of transcripts from which I want to extract sequences
cut -f4 ../background_transcriptome/mm10_GRCm38_MM257.tpm1_coding.utr3.exons_most.expressed.coding.isoform.per.gene.bed | awk 'BEGIN {FS="_"} {print $1"_"$2}' | sort | uniq >| temp.Refseq.ids
bsub -q 18 -K Rscript ${RDIR}/obtain_3UTRseq_from_bed_190220.R utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line.bed temp.Refseq.ids 17 utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line sleep 10 &
wait


### remove line breaks in fasta 
awk '/^>/{print s? s"\n"$0:$0;s="";next}{s=s sprintf("%s",$0)}END{if(s)print s}' utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line.fa >| temp
mv -f temp utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line.fa 

### generate randomly selected background sequences of 35bp from DAZL bound 3' UTRs that are express at a min TPM of 1; most robustly expressed isoform per gene
### get motif frequence relative to center position in background sequences
### background comprise of randomly selected 35nt sequences from 3' UTR of DAZL bound transcript (canonical transcripts only, as identified by UCSC)
mkdir background_DAZL_bound_random_35bp
cd background_DAZL_bound_random_35bp

ln -s ../../background_transcriptome/3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound.fa  . 

### randomly generate 20 35bp sequences per 3' UTR that is DAZL bound
python ${PYTHDIR}/random_selection_35bp_seq.py 3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound.fa random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound.fa


### get a stats file and graph for each triplicate motif
cd ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/metaplot_around_CL_nt
for FIRST in U G C A      ### Outer for loop ###
do

    for SECOND in U G C A ### Middle inner for loop ###
    do
         
        for THIRD in U G C A ### Inner inner for loop ###
        do
            MOTIF=${FIRST}${SECOND}${THIRD}
            echo ${MOTIF}
            ### identify all instances of a specific motif in the randomly generated sequences (i.e., background file)
            BKGDDIR="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/20180730_ASPeak/ASPeak_IgG_RNAseq_gapno0_CL.site.only_w.intergenic/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/metaplot_around_CL_nt/background_DAZL_bound_random_35bp"
            cd ${BKGDDIR}
            python ${PYTHDIR}/identify_all_tri.nt_motifs_fa_input_only_180409.py ${MOTIF} ${BKGDDIR}/random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound.fa random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound_${MOTIF}.bed
            ### get motif frequency relative to each crosslinked site
            cd ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/metaplot_around_CL_nt
            python ${PYTHDIR}/identify_all_tri.nt_motifs_fa_input_only_180409.py ${MOTIF} utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line.fa utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line_${MOTIF}.bed
            cut -f6 utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line_${MOTIF}.bed | sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' >| CL_${MOTIF}_tmp
            cut -f6 background_DAZL_bound_random_35bp/random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound_${MOTIF}.bed| sort | uniq -c  | sed -r 's/^( *[^ ]+) +/\1\t/'  >| bkgd_${MOTIF}_tmp
            CLseq=$(wc -l < utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line.bed) ### capture number of CL sites in variale CLseq
            temp=$(wc -l < background_DAZL_bound_random_35bp/random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound.fa) #### capture number of background sites in variable bkgdseq
            bkgdseq=$(expr $temp / 2)
            join -1 2 -2 2 -t $'\t' CL_${MOTIF}_tmp bkgd_${MOTIF}_tmp | awk -v a="$CLseq" -v b="$bkgdseq"  'BEGIN {FS="\t"} {print $1"\t"$2"\t"$3"\t"a"\t"b}' | awk 'BEGIN{print "nt.from.CL.site\tno.motifs.CL\tno.motifs.bkgd\tno.CL.sites\tno.bkgd.seq"}1' >| CL_bkgd_${MOTIF}_counts
            rm -f *tmp
            Rscript ${RDIR}/crosslinked_metaplots.R ${MOTIF}
        done

    done

done



### get a stats file and graph for additional motifs
cd ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/metaplot_around_CL_nt
for MOTIF in UGUU UGUU.notG UGUU.UorA UUU.CorG.UUU GUUG GUUC; do 
    ### identify all instances of a specific motif in the randomly generated sequences (i.e., background file)
    BKGDDIR="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/20180730_ASPeak/ASPeak_IgG_RNAseq_gapno0_CL.site.only_w.intergenic/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/metaplot_around_CL_nt/background_DAZL_bound_random_35bp"
    cd ${BKGDDIR}
    python ${PYTHDIR}/identify_all_${MOTIF}_fa_input_only_180409.py ${BKGDDIR}/random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound.fa random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound_${MOTIF}.bed
    ### get motif frequence relative to each crosslinked site
    cd ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/metaplot_around_CL_nt
    python ${PYTHDIR}/identify_all_${MOTIF}_fa_input_only_180409.py utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line.fa utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line_${MOTIF}.bed
    cut -f6 utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line_${MOTIF}.bed | sort | uniq -c | sed -r 's/^( *[^ ]+) +/\1\t/' >| CL_${MOTIF}_tmp
    cut -f6 background_DAZL_bound_random_35bp/random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound_${MOTIF}.bed | sort | uniq -c  | sed -r 's/^( *[^ ]+) +/\1\t/'  >| bkgd_${MOTIF}_tmp
    CLseq=$(wc -l < utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line.bed) ### capture number of CL sites in variale CLseq
    temp=$(wc -l < background_DAZL_bound_random_35bp/random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound.fa) #### capture number of background sites in variable bkgdseq
    bkgdseq=$(expr $temp / 2)
    join -1 2 -2 2 -t $'\t' CL_${MOTIF}_tmp bkgd_${MOTIF}_tmp | awk -v a="$CLseq" -v b="$bkgdseq"  'BEGIN {FS="\t"} {print $1"\t"$2"\t"$3"\t"a"\t"b}' | awk 'BEGIN{print "nt.from.CL.site\tno.motifs.CL\tno.motifs.bkgd\tno.CL.sites\tno.bkgd.seq"}1' >| CL_bkgd_${MOTIF}_counts
    rm -f *tmp
    Rscript ${RDIR}/crosslinked_metaplots.R ${MOTIF}
done



### plot for publication: one graph with motifs related to my major motif: GUU, UUU, UGUU(U/A)
### another graph with motifs related to my major motif: UGUU, GUU(U/A), UGUU(U/A)
## one graph for motifs that were previously reported: UGUU(U/A) vs. GUUG (Zagore et al., 2018), GUUC (Maegawa et al., 2002 Genes to Cells; Reynolds et al., 20015 Human Mol Genetics), UUU[C/G]UUU) (Chen et al 2011 Genes and Dev Conti lab)
bsub Rscript ${RDIR}/crosslinked_metaplots_for_publication.R 


###kpLogo analysis
cd ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/
mkdir kpLogo
cd kpLogo
ln -s ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/metaplot_around_CL_nt/utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line.bed .

	### using bed file that contains crosslinked position as a single entry (i.e., consecutive crosslinked positions are in separate lines)
	###reorganize and clean up bed file; shorten name so it only corresponds to one transcript per crosslinked site
awk 'BEGIN {FS="\t"} {print $1"\t"$2"\t"$3"\t"$12"\t"$6"\t"$4} ' utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_one.CL.site.per.line.bed | awk 'BEGIN {FS=","} {print $1} ' | awk 'BEGIN {FS="\t"} {print $1"\t"$2"\t"$3"\t"$6"\t"$4"\t"$5}' >| temp.bed


### expand bed file by 10 nucleotides on each side of crosslinked sites; extract as fasta sequences
	### "-q 18" send job to ubuntu cluster 18 nodes
### R script requries a list of transcripts from which I want to extract sequences
cut -f4 ../background_transcriptome/mm10_GRCm38_MM257.tpm1_coding.utr3.exons_most.expressed.coding.isoform.per.gene.bed | awk 'BEGIN {FS="_"} {print $1"_"$2}' | sort | uniq >| temp.Refseq.ids
bsub -q 18 -K Rscript ${RDIR}/obtain_3UTRseq_from_bed_kpLogo.R temp.bed temp.Refseq.ids 10 temp_CLsites.10bp sleep 10 &
wait

	### get closest GUU for each crosslinked site
	### remove sequences that are shorter than 21nt (because crosslinked site is too close to the beginning or end of 3' UTR)
python ${PYTHDIR}/identify_closest_5GUU3_190218.py temp_CLsites.10bp.fa temp.bed utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_w.closest.GUU


	### create bed file containing only closest GUU per crosslinked site; remove "n/a"s (crosslinked sites with no GUU found; remove repeated lines
	### file output with -log10(pval) in name position
awk '{if ($7 != "n/a") print $1"\t"$6"\t"$7"\t"$10"\tn/a\t"$4}' utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_w.closest.GUU | sort -k1,1 -k2,2n | uniq >| GUU_closest_to_CL_sites.bed

###get fasta file of closest GUU to crosslinked site, plus 10nt and minus 10 nt (23nt total)
bsub -q 18 -K Rscript ${RDIR}/obtain_3UTRseq_from_bed_190220.R GUU_closest_to_CL_sites.bed temp.Refseq.ids 10 GUU_closest_to_CL_sites sleep 10 &
wait
rm -f temp*
	

### reorganize fasta file so that it is in the following format: sequence -log10(pval)
sed 'N;s/\n/ /' GUU_closest_to_CL_sites.fa | awk 'BEGIN {FS=" ";} {print $1"::"$2}' |  awk 'BEGIN {FS="::";} {print $2"\t"$1}' |awk 'BEGIN {FS=">";} {print $1$2}' >| kpLogo_input

### obtain background sequences for kpLogo analysis
### background sequences are all GUUs in the 3'UTRs of Dazl targets
mkdir background_DAZL_bound_utr3
cd background_DAZL_bound_utr3
ln -s ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/metaplot_around_CL_nt/background_DAZL_bound_random_35bp/random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound.fa .
	### extract GUUs (plus and minus 10nt of either side of motif; 23nt in total) from background sequences
python ${PYTHDIR}/identify_all_GUU_23nt_output_190204.py random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound.fa random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound_GUUs.txt
cut -f2 random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound_GUUs.txt | sort  >| temp_background
cut -f2 random_35bp_seqs_from_3UTRs_most.expressed.coding.isoform.per.gene_DAZL.bound_GUUs.txt | sort | awk 'BEGIN {FS="\t"} {print $1"\t0"}' >| temp_background


### add background sequences to target list
cd ../
cat background_DAZL_bound_utr3/temp_background >> kpLogo_input
sort -k1,1n kpLogo_input >| kpLogo_input_sorted


### running kpLogo with weighted setting centered around GUU
### for the GUUs closest to crosslinked sites, the score is the -log10(pval) at crosslinked site
### for background sequence, score of 0
~/bin/kpLogo/bin/kpLogo kpLogo_input_sorted -weighted 





######################################################################
#####
### Determine whether DAZL binding sites are conserved among vertebrates
#####

cd ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis

mkdir conservation_analysis 
cd conservation_analysis

# bed file containing the closest GUU to crosslinked sites
ln -s ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/kpLogo/GUU_closest_to_CL_sites.bed . 

# sort bed file of Dazl binding sites in 3' UTR; removed extra columns
bedtools sort -i ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1.bed | awk 'BEGIN {FS="\t"} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6}' >|  utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_sorted.bed

# Link Dazl bound transcripts' 3' UTRs (longest expressed isoform per transcript)
ln -s ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/background_transcriptome/mm10_GRCm38_MM257.tpm1_coding.utr3.exons_most.expressed.coding.isoform.per.gene_DAZL.bound.bed .

# Create a background file that contains Dazl bound 3'UTRs without DAZL-bound nucleotides
bedtools subtract -s -a mm10_GRCm38_MM257.tpm1_coding.utr3.exons_most.expressed.coding.isoform.per.gene_DAZL.bound.bed -b utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_sorted.bed >| mm10_GRCm38_MM257.tpm1_coding.utr3.exons_most.expressed.coding.isoform.per.gene_DAZL.bound_no.Dazl.bound.nts.bed

# Download Vertebrate files for PhyloP scores and PhastCon scores from UCSC Genome Browser
# Convert bw file to bedGraph file
# Intersect downloaded bedGraph file wtih DAZL binding sites in 3' UTR


mkdir phastCons 
cd phastCons
bsub "rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/README.txt ."
bsub -K "rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/phastCons60way/mm10.60way.phastCons.bw ."
wait 
bsub -K "bigWigToBedGraph mm10.60way.phastCons.bw mm10.60way.phastCons.bedGraph"
wait
bsub -K "gzip -f mm10.60way.phastCons.bw"
wait
rm -f mm10.60way.phastCons.bw
#cat mm10.60way.phastCons.bedGraph | awk 'BEGIN {FS="\t"} {print $1"\t"$2"\t"$3"\tname\t"$4"\t."}' >| mm10.60way.phastCons.bed
bsub "bedtools intersect -a mm10.60way.phastCons.bedGraph -b ../utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_sorted.bed >| ../utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_phastCons.bed"
#bsub "bedtools intersect -a mm10.60way.phastCons.bedGraph -b ../GUU_closest_to_CL_sites.bed >| ../GUU_closest_to_CL_sites_phastCons.bed"
bsub "bedtools intersect -a mm10.60way.phastCons.bedGraph -b ../mm10_GRCm38_MM257.tpm1_coding.utr3.exons_most.expressed.coding.isoform.per.gene_DAZL.bound_no.Dazl.bound.nts.bed >| ../mm10_GRCm38_MM257.tpm1_coding.utr3.exons_most.expressed.coding.isoform.per.gene_DAZL.bound_no.Dazl.bound.nts_phastCons.bed"


cd /lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/20180730_ASPeak/ASPeak_IgG_RNAseq_gapno0_CL.site.only_w.intergenic/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/conservation_analysis
mkdir phyloP
cd phyloP
bsub -K "rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/phyloP60way/README.txt ."
bsub -K "rsync -avz --progress rsync://hgdownload.cse.ucsc.edu/goldenPath/mm10/phyloP60way/mm10.60way.phyloP60way.bw ."
wait
bsub -K "bigWigToBedGraph mm10.60way.phyloP60way.bw mm10.60way.phyloP60way.bedGraph"
wait
bsub -K "gzip -f mm10.60way.phyloP60way.bw"
wait
rm -f mm10.60way.phyloP60way.bw
#cat mm10.60way.phyloP60way.bedGraph | awk 'BEGIN {FS="\t"} {print $1"\t"$2"\t"$3"\tname\t"$4"\t."}' >| mm10.60way.phyloP60way.bed
bsub "bedtools intersect -a mm10.60way.phyloP60way.bedGraph -b ../utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_sorted.bed >| ../utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_phyloP.bed"
#bsub "bedtools intersect -a mm10.60way.phyloP60way.bedGraph -b ../GUU_closest_to_CL_sites.bed >| ../GUU_closest_to_CL_sites_phyloP60way.bed"
bsub -K "bedtools intersect -a mm10.60way.phyloP60way.bedGraph -b ../mm10_GRCm38_MM257.tpm1_coding.utr3.exons_most.expressed.coding.isoform.per.gene_DAZL.bound_no.Dazl.bound.nts.bed >| ../mm10_GRCm38_MM257.tpm1_coding.utr3.exons_most.expressed.coding.isoform.per.gene_DAZL.bound_no.Dazl.bound.nts_phyloP.bed"
wait

## will run Rscript code at end to get graphs and stats for these datasets



#############
### compare the conservation of the most robust DAZL crosslinked nucleotides in each gene to other DAZL crosslinked nucleotides
#############
cd ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/conservation_analysis

mkdir strongest_peaks
cd strongest_peaks 

# link bed file of DAZL binding sites
ln -s  ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1.bed 


# obtain a bed file containing the most robust peak (based on pval) per gene; if there are multiple peaks with equal pvalues tied for most robust peak, pick the one with the most reads from CLIP data
sort -k4,4 -k12,12n -k5,5nr utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1.bed | awk 'BEGIN {FS="\t"} {print $0"\t"$4}' | uniq -f13 | awk 'BEGIN {FS="\t"} {print $1"\t"$2"\t"$3"\t"$4"\t"$5"\t"$6"\t"$7"\t"$8"\t"$9"\t"$10"\t"$11"\t"$12"\t"$13}' | bedtools sort >| utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_most.robust.peak.per.gene.bed

grep -wvFf utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_most.robust.peak.per.gene.bed utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1.bed | bedtools sort >| utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_not.most.robust.peak.per.gene.bed


for FILE in utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_most.robust.peak.per.gene utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_not.most.robust.peak.per.gene; do 
bsub "bedtools intersect -a ../utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_phastCons.bed -b ${FILE}.bed >| ${FILE}_phastCons.bed"
bsub "bedtools intersect -a ../utr3_coding.only.peaks_present.in.at.least.2.replicates_min.TPM.1_phyloP.bed  -b ${FILE}.bed >| ${FILE}_phyloP.bed"
done

##graphs and stats for CL vs noncrosslinked analysis; also for most robust CL sites vs other CL sites
cd ${DIR}/foundpeaks/merged.peaks.190622/utr3_coding.only/min.TPM.1/sequence_analysis/conservation_analysis
bsub Rscript ${RDIR}/conservation_graphs.R









#### Identify whether GUU motif is enriched in DAZL binding sites for each region
for REGION in cds ncRNA.no.rRNA retrogenes utr5_coding.only utr3_coding.only; do
	cd ${DIR}/foundpeaks/merged.peaks.190622/${REGION}/min.TPM.1/sequence_analysis/
	mkdir GUU_analysis_ame
	cd GUU_analysis_ame
	ln -s ../../../${REGION}.peaks_present.in.at.least.2.replicates.bed


	### expand bed file by 2 nucleotides on each side of crosslinked sites for GUU analysis; add line number to each name so each sequences has a uniqueID using NR
	awk '{FS="\t"} {if (NR==1) {print $0}  else {print $1"\t"$2-2"\t"$3+2"\t"$4"_"NR "\t"$5"\t"$6}} ' ${REGION}.peaks_present.in.at.least.2.replicates.bed >| temp_expanded_GUU.bed 

	### expand bed file by 3 nucleotides on each side of crosslinked sites for UGUU(UA) analysis; add line number to each name so each sequences has a uniqueID using NR
	awk '{FS="\t"} {if (NR==1) {print $0}  else {print $1"\t"$2-3"\t"$3+3"\t"$4"_"NR "\t"$5"\t"$6}} ' ${REGION}.peaks_present.in.at.least.2.replicates.bed >| temp_expanded_UGUU.UorA.bed 


	### get fasta file
	FA="/nfs/genomes/mouse_mm10_dec_11_no_random/fasta_whole_genome/mm10.fa"
	bedtools getfasta -name -s -fi ${FA} -bed temp_expanded_GUU.bed  >| temp_expanded_GUU.fa
	bedtools getfasta -name -s -fi ${FA} -bed temp_expanded_UGUU.UorA.bed  >| temp_expanded_UGUU.UorA.fa

	## replace T with U
	sed 's/T/U/g'  temp_expanded_GUU.fa>| temp_expanded_GUU_rna.seq.fa
	sed 's/T/U/g'  temp_expanded_UGUU.UorA.fa>| temp_expanded_UGUU.UorA_rna.seq.fa


	### shorten header names so that meme can deal with them more easily
	awk 'BEGIN {FS="\n"} {if (substr($1,1,1) == ">") print ">"NR; else print $0;} ' temp_expanded_GUU_rna.seq.fa>| temp_expanded_GUU_rna.seq_header.modified.fa
	awk 'BEGIN {FS="\n"} {if (substr($1,1,1) == ">") print ">"NR; else print $0;} ' temp_expanded_UGUU.UorA_rna.seq.fa>| temp_expanded_UGUU.UorA_rna.seq_header.modified.fa


	### using RNA library using pre-created RNA alphabet file in meme format http://meme-suite.org/doc/alphabet-format.html?man_type=web
	yes| cp  -i ~/bin/for_meme_motif_analysis/alphabet_rna.txt alphabet_rna.txt 
	### get GUU, UGUU(U/A) motif in meme format
	### using RNA library using pre-created RNA alphabet file in meme format http://meme-suite.org/doc/alphabet-format.html?man_type=web
	/usr/local/meme/bin/iupac2meme -alph alphabet_rna.txt GUU >| GUU_motif_meme_format
	/usr/local/meme/bin/iupac2meme -alph alphabet_rna.txt UGUUW >| UGUU.UorA_motif_meme_format
	
	## file for both motifs at once
	/usr/local/meme/bin/iupac2meme -alph alphabet_rna.txt GUU UGUUW >| motifs_motif_meme_format



	### create background of shuffled sequences
	/usr/local/meme/bin/fasta-shuffle-letters temp_expanded_GUU_rna.seq_header.modified.fa shuffled_control_GUU.fa
	/usr/local/meme/bin/fasta-shuffle-letters temp_expanded_UGUU.UorA_rna.seq_header.modified.fa shuffled_control_UGUU.UorA.fa

	### search for GUU, UGUU(U/A) motif using ame; set max p value to 1 so I can see non significant enrichment
	/usr/local/meme/bin/ame --verbose 1 --oc . --control shuffled_control_GUU.fa --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 1 temp_expanded_GUU_rna.seq_header.modified.fa GUU_motif_meme_format
	mv -f ame.html ${REGION}_GUU_ame.html
	/usr/local/meme/bin/ame --verbose 1 --oc . --control shuffled_control_UGUU.UorA.fa --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 1 temp_expanded_UGUU.UorA_rna.seq_header.modified.fa UGUU.UorA_motif_meme_format
	mv -f ame.html ${REGION}_UGUU.UorA_ame.html
	
	## analyze both motifs at once BUT using longer input sequence (not ideal for shorter GUU motif)
	/usr/local/meme/bin/ame --verbose 1 --oc . --control shuffled_control_UGUU.UorA.fa --bgformat 1 --scoring avg --method ranksum --pvalue-report-threshold 1 temp_expanded_UGUU.UorA_rna.seq_header.modified.fa motifs_motif_meme_format
	mv -f ame.html ${REGION}_motifs_ame.html

	##make motif for GUU, UGUU(U/A)
	/usr/local/meme/bin/ceqlogo -i1 GUU_motif_meme_format -o GUU.png -f PNG
	/usr/local/meme/bin/ceqlogo -i1 UGUU.UorA_motif_meme_format -o UGUU.UorA.png -f PNG
	
	#rm -f temp*

done
