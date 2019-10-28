# !/bin/bash

### Maria Mikedis mikedis@wi.mit.edu


cd /lab/solexa_page/maria/
mkdir iCLIP_DAZL_undiff.gonia
cd iCLIP_DAZL_undiff.gonia

mkdir 170314_iCLIP_rawdata
cd 170314_iCLIP_rawdata

cp -r  /lab/solexa_public/Page/170309_WIGTC-HISEQ2A_CAPK4ANXX .
cd QualityScore/
gunzip -c unknown-s_6_1_sequence.txt.tar.gz | tar xopf - 
gzip unknown-s_6_1_sequence.txt.gz


### quality cutting; remove any 3' adapter that may be in read

cd ../..
mkdir 170314_quality_trimmed
cd 170314_quality_trimmed
ln -s /lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/170314_iCLIP_rawdata/170309_WIGTC-HISEQ2A_CAPK4ANXX/QualityScore/unknown-s_6_1_sequence.txt.gz .

gzip -d /lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/170314_iCLIP_rawdata/170309_WIGTC-HISEQ2A_CAPK4ANXX/QualityScore/unknown-s_6_1_sequence.txt.gz 


bsub "fastq_quality_filter -i /lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/170314_iCLIP_rawdata/170309_WIGTC-HISEQ2A_CAPK4ANXX/QualityScore/unknown-s_6_1_sequence.txt -Q 64  -v -q 20 -p 75 -o unknown-s_6_1_sequence.trimmed.txt"

mv unknown-s_6_1_sequence.trimmed.txt.gz unknown-s_6_1_sequence.trimmed.fastq.gz

cutadapt -a AGATCGGAAGAGCGGTTCAG -m 24 -o unknown-s_6_1_sequence.trimmed2.txt unknown-s_6_1_sequence.trimmed.fastq.gz

mv unknown-s_6_1_sequence.trimmed.fastq.gz unknown-s_6_1_sequence.trimmed.txt.gz

gzip /lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/170314_iCLIP_rawdata/170309_WIGTC-HISEQ2A_CAPK4ANXX/QualityScore/unknown-s_6_1_sequence.txt





### remove PCR duplicates

cd ..
mkdir 170314_removal_PCR_duplicates
cd 170314_removal_PCR_duplicates


bsub "fastx_collapser -v -i  /lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/170314_quality_trimmed/unknown-s_6_1_sequence.trimmed2.txt > unknown-s_6_1_sequence.trimmed2.PCR.dup.removed.fasta"

bsub "gzip /lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/170314_quality_trimmed/unknown-s_6_1_sequence.trimmed2.txt"





### demultiplex
cd ..
mkdir 170314_demultiplex
cd 170314_demultiplex

INPUT="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/170314_removal_PCR_duplicates/unknown-s_6_1_sequence.trimmed2.PCR.dup.removed.fasta"

	###remove first three nucleotides (random barcode); trim 3' adapter, if any is on read

fastx_trimmer -f 4 -i ${INPUT} -o unknown-s_6_1_sequence.trimmed2.PCR.dup.removed_ran.bc.trim.fasta

cat unknown-s_6_1_sequence.trimmed2.PCR.dup.removed_ran.bc.trim.fasta | fastx_barcode_splitter.pl --bcfile barcodes --prefix demultiplexed_ --bol 

gzip unknown-s_6_1_sequence.trimmed2.PCR.dup.removed_ran.bc.trim.fasta

	### trim off regular barcodes (4 nucleotides) and 2 nucleotides of random barcode on  5' end

for FILE in demultiplexed_MM218_2_IgG_CLIP demultiplexed_MM222_2_2_IgG_CLIP demultiplexed_MM218_3_DAZL_noCL demultiplexed_MM222_2_3_DAZL_noCL demultiplexed_MM218_4_input demultiplexed_MM222_2_4_input demultiplexed_MM222_1_1_DAZL_CLIP demultiplexed_MM222_1_2_IgG_CLIP demultiplexed_MM222_1_3_DAZL_noCL; do 


	fastx_trimmer -f 7 -i ${FILE} -o ${FILE}_trimmed

done


cd /lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia
mkdir 180411_finalized_analysis


DIR="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis"

### run STAR for iCLIP reads wtih 30nt overhang
### modify output bam files so only crosslinked nt is represented by each read
### index bam files

cd ${DIR}
now1="$(date +'%Y%m%d')"
mkdir ${now1}_iCLIP_STAR
CLIPDIR="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/${now1}_iCLIP_STAR"
cd ${CLIPDIR}

ln -s /lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/170314_demultiplex/*trimmed .


INDEX="/lab/solexa_page/maria/genome_files/Mus.musculus_GRCm38.p5_rel90/STAR/STAR_output_overhang30/"
jobnum=1

for FILE in demultiplexed_MM218_1_DAZL_CLIP_trimmed demultiplexed_MM218_2_IgG_CLIP_trimmed demultiplexed_MM218_3_DAZL_noCL_trimmed demultiplexed_MM218_4_input_trimmed demultiplexed_MM222_1_1_DAZL_CLIP_trimmed demultiplexed_MM222_1_2_IgG_CLIP_trimmed demultiplexed_MM222_1_3_DAZL_noCL_trimmed demultiplexed_MM222_1_4_input_trimmed demultiplexed_MM222_2_1_DAZL_CLIP_trimmed demultiplexed_MM222_2_2_IgG_CLIP_trimmed demultiplexed_MM222_2_3_DAZL_noCL_trimmed demultiplexed_MM222_2_4_input_trimmed; do 

	bsub -J "iCLIP.jobs${jobnum}" -n8 -R "span[hosts=1]" "/nfs/apps/STAR/2.5.4b/STAR --version STAR_2.5.4b --genomeDir ${INDEX} --readFilesIn ${FILE}  --outFileNamePrefix ${FILE}_STAR_ --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --outFilterMismatchNmax 2 --runThreadN 8 --outSAMattributes None --outReadsUnmapped Fastx"

	bsub -w "done(iCLIP.jobs${jobnum})" -J "iCLIP.jobs${jobnum}.1" "python /lab/solexa_page/maria/scripts/CLIP_analysis/converting_mapped_reads_to_crosslinked_sites.py ${FILE}_STAR_Aligned.out.sam ${FILE}_STAR_Aligned.out_CL.site.only.sam"

	bsub -w "done(iCLIP.jobs${jobnum}.1)" -J "iCLIP.jobs${jobnum}.2" "samtools view -bT /lab/solexa_page/maria/genome_files/Mus.musculus_GRCm38.p5_rel90/mm10.fa  ${FILE}_STAR_Aligned.out.sam >|  ${FILE}_STAR_Aligned.out.bam"
	bsub -w "done(iCLIP.jobs${jobnum}.1)" -J "iCLIP.jobs${jobnum}.3" "samtools view -bT /lab/solexa_page/maria/genome_files/Mus.musculus_GRCm38.p5_rel90/mm10.fa  ${FILE}_STAR_Aligned.out_CL.site.only.sam >|  ${FILE}_STAR_Aligned.out_CL.site.only.bam"

	bsub -w "done(iCLIP.jobs${jobnum}.2)" -J "iCLIP.jobs${jobnum}.4" "rm -f ${FILE}_STAR_Aligned.out.sam"
	bsub -w "done(iCLIP.jobs${jobnum}.3)" -J "iCLIP.jobs${jobnum}.5" "rm -f ${FILE}_STAR_Aligned.out_CL.site.only.sam"

	bsub  -w "done(iCLIP.jobs${jobnum}.4)" -J "iCLIP.jobs${jobnum}.6" "samtools sort ${FILE}_STAR_Aligned.out.bam >| ${FILE}_STAR_Aligned.out_sorted.bam"
	bsub  -w "done(iCLIP.jobs${jobnum}.5)" -J "iCLIP.jobs${jobnum}.7" "samtools sort ${FILE}_STAR_Aligned.out_CL.site.only.bam >| ${FILE}_STAR_Aligned.out_CL.site.only_sorted.bam" 

	bsub  -w "done(iCLIP.jobs${jobnum}.6)" -J "iCLIP.jobs${jobnum}.8" "samtools index ${FILE}_STAR_Aligned.out_sorted.bam ${FILE}_STAR_Aligned.out_sorted.bam.bai"
	bsub  -w "done(iCLIP.jobs${jobnum}.7)" -J "iCLIP.jobs${jobnum}.9" "samtools index ${FILE}_STAR_Aligned.out_CL.site.only_sorted.bam ${FILE}_STAR_Aligned.out_CL.site.only_sorted.bam.bai"

	jobnum=$((jobnum + 1))

done

### run STAR for 3S RNAseq reads wtih 39nt overhang
cd ${DIR}
now2="$(date +'%Y%m%d')"
mkdir ${now2}_3S_RNAseq_STAR
RNADIR="/lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/180411_finalized_analysis/${now2}_3S_RNAseq_STAR"
cd ${now2}_3S_RNAseq_STAR


ln -s /lab/solexa_page/maria/iCLIP_DAZL_undiff.gonia/171014_undiff_gonia_3S_rnaseq/171014_rawdata/170928_WIGTC-HISEQ2B_CB8HMANXX/QualityScore/lab/solexa_public/Page/170928_WIGTC-HISEQ2B_CB8HMANXX/QualityScore/*.gz . .


INDEX="/lab/solexa_page/maria/genome_files/Mus.musculus_GRCm38.p5_rel90/STAR/STAR_output_overhang39/"
jobnum=1
for FILE in CTTGTA-s_2_1_sequence GGCTAC-s_2_1_sequence; do 
	
	bsub -J "RNA.jobs${jobnum}" -n8 -R "span[hosts=1]" "STAR --genomeDir ${INDEX} --readFilesIn ${FILE}.txt.gz --readFilesCommand zcat --outFileNamePrefix ${FILE}_STAR_ --outFilterMultimapNmax 1 --alignEndsType Extend5pOfRead1 --outFilterMismatchNmax 2 --runThreadN 8 --outReadsUnmapped Fastx"

	bsub -w "done(RNA.jobs${jobnum})" -J "RNA.jobs${jobnum}.1" "samtools view -bT /lab/solexa_page/maria/genome_files/Mus.musculus_GRCm38.p5_rel90/mm10.fa ${FILE}_STAR_Aligned.out.sam >|  ${FILE}_STAR_Aligned.out.bam"

	bsub -w "done(RNA.jobs${jobnum}.1)" -J "RNA.jobs${jobnum}.2" "rm -f ${FILE}_STAR_Aligned.out.sam"

	bsub  -w "done(RNA.jobs${jobnum}.1)" -J "RNA.jobs${jobnum}.3" "samtools sort ${FILE}_STAR_Aligned.out.bam >| ${FILE}_STAR_Aligned.out_sorted.bam"

	bsub  -w "done(RNA.jobs${jobnum}.3)" -J "RNA.jobs${jobnum}.4" "samtools index ${FILE}_STAR_Aligned.out_sorted.bam ${FILE}_STAR_Aligned.out_sorted.bam.bai"

	jobnum=$((jobnum + 1))

done




### ASPeak analysis, using crosslinked nucleotides only for CLIP files; IgG CLIP as background; ASPeak input bed files that include intergenic regions


cd ${DIR}
now3="$(date +'%Y%m%d')"
mkdir ${now3}_ASPeak
cd ${now3}_ASPeak

### link files for mapped RNAseq reads; will be broken until preceding jobs are done
for FILE in CTTGTA-s_2_1_sequence GGCTAC-s_2_1_sequence; do 
bsub -w "done(RNA.jobs*)"  "ln -s ${RNADIR}/${FILE}_STAR_Aligned.out_sorted.bam ."
done

### link files for mapped CLIP reads; will be broken until preceding jobs are done
for FILE in demultiplexed_MM218_1_DAZL_CLIP_trimmed demultiplexed_MM218_2_IgG_CLIP_trimmed demultiplexed_MM218_3_DAZL_noCL_trimmed demultiplexed_MM218_4_input_trimmed demultiplexed_MM222_1_1_DAZL_CLIP_trimmed demultiplexed_MM222_1_2_IgG_CLIP_trimmed demultiplexed_MM222_1_3_DAZL_noCL_trimmed demultiplexed_MM222_1_4_input_trimmed demultiplexed_MM222_2_1_DAZL_CLIP_trimmed demultiplexed_MM222_2_2_IgG_CLIP_trimmed demultiplexed_MM222_2_3_DAZL_noCL_trimmed demultiplexed_MM222_2_4_input_trimmed; do 
bsub -w "done(iCLIP.jobs*)" "ln -s ${CLIPDIR}/${FILE}_STAR_Aligned.out_CL.site.only_sorted.bam ."
done


## call peaks using crosslinked nucleotide only
CLIP1=demultiplexed_MM218_1_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted.bam
CLIP2=demultiplexed_MM222_1_1_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted.bam
CLIP3=demultiplexed_MM222_2_1_DAZL_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted.bam
IgG1=demultiplexed_MM218_2_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted.bam
IgG2=demultiplexed_MM222_1_2_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted.bam
IgG3=demultiplexed_MM222_2_2_IgG_CLIP_trimmed_STAR_Aligned.out_CL.site.only_sorted.bam
INPUT1=CTTGTA-s_2_1_sequence_STAR_Aligned.out_sorted.bam
INPUT2=GGCTAC-s_2_1_sequence_STAR_Aligned.out_sorted.bam


### rerunning analysis after adding intergenic bed file to ASPeak files on 3/26/18
bsub -w "done(iCLIP.jobs*) && done(RNA.jobs*)" -J "ASPeak" "ASPeak.pl -param ~/bin/ASPeak/scripts/default_parameters_mm10.txt -gapnumber 0 -bed /lab/solexa_page/maria/genome_files/ASPeak_171122/bed_files_modified_for_ASPeak -lib ${CLIP1}:${CLIP2}:${CLIP3} -rnaseq ${INPUT1}:${INPUT2} -outdir ASPeak_IgG_RNAseq_gapno0_CL.site.only_w.intergenic -control ${IgG1}:${IgG2}:${IgG3}"
