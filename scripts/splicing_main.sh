#!/bin/bash

outDir=$1
pData=$2
gtf=$3
condPairs=$4
nthreads=$5
indexAppendix=$6
readLength=$7

dexseq_script="/home/scripts/dexAna.R"
rmats_script="/home/software/rMATS.4.0.2/rMATS-turbo-Linux-UCS4/rmats.py"
drimseq_script="/home/scripts/run_DRIMSeq.R"

mkdir -p $outDir/diff_splicing_outs/rMATS/

#DEXSeq
( [ -f "$outDir/diff_splicing_outs/DEXSeq.out" ] && echo "[INFO] [DEXSeq] $outDir/diff_splicing_outs/DEXSeq.out already exists, skipping.."$'\n' ) || ( echo "[INFO] [DEXSeq] ["`date "+%Y/%m/%d-%H:%M:%S"`"] $dexseq_script $pData $outDir/COUNTS/DEXSeq_HTcounts $condPairs /home/indices/DEXSeq/$indexAppendix/annot.noaggregate.gtf $outDir/diff_splicing_outs" && $dexseq_script $pData $outDir/COUNTS/DEXSeq_HTcounts $condPairs /home/indices/DEXSeq/$indexAppendix/annot.noaggregate.gtf $outDir/diff_splicing_outs && echo "[INFO] [DEXSeq] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Splicing analysis finished"$'\n' )



#rMATS
file1=$outDir/diff_splicing_outs/rMATS/b1.txt
file2=$outDir/diff_splicing_outs/rMATS/b2.txt
#[ ! -f $file1 ] && [ ! -f $file2 ] && touch $file1 $file2
#while read line; do 
#	cond=`cut -f2 line`
#	[ "$cond" == "0" ] && 
#	[ "$cond" == "1" ] && 
#done < $pData

echo "[INFO] [rMATS] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Creating b-files..."
cat $pData | grep -P "\t0" | awk '$0="/home/data/out/HISAT/dta/"$0'  | cut -f1 |  tr "\\n" "," | sed 's/,/.bam,/'g | rev | cut -c 2- | rev > $file1
cat $pData | grep -P "\t1" | awk '$0="/home/data/out/HISAT/dta/"$0'  | cut -f1 |  tr "\\n" "," | sed 's/,/.bam,/'g | rev | cut -c 2- | rev > $file2

( [ -f "$outDir/diff_splicing_outs/rMATS/SE.MATS.JC.txt" ] && echo "[INFO] [rMATS] $outDir/diff_splicing_outs/rMATS/SE.MATS.JC.txt already exists, skipping.."$'\n' ) || ( echo "[INFO] [rMATS] ["`date "+%Y/%m/%d-%H:%M:%S"`"] python $rmats_script --b1 $outDir/diff_splicing_outs/rMATS/b1.txt --b2 $outDir/diff_splicing_outs/rMATS/b2.txt -t paired --libType fr-unstranded --readLength $readLength --gtf $gtf --od $outDir/diff_splicing_outs/rMATS --nthread $nthreads" && python $rmats_script --b1 $outDir/diff_splicing_outs/rMATS/b1.txt --b2 $outDir/diff_splicing_outs/rMATS/b2.txt\
 	-t paired --libType fr-unstranded --readLength $readLength --gtf $gtf --od $outDir/diff_splicing_outs/rMATS --nthread $nthreads )

echo "[INFO] [rMATS] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Splicing analysis finished"$'\n'



#DRIMSeq
( [ -f "$outDir/diff_splicing_outs/DRIMSeq.out" ] && echo "[INFO] [DRIMSeq] $outDir/diff_splicing_outs/DRIMSeq.out already exists, skipping.."$'\n' ) || ( echo "[INFO] [DRIMSeq] ["`date "+%Y/%m/%d-%H:%M:%S"`"] $drimseq_script --counts $outDir/KALLISTO/alignment --pdata $pData --outfile $outDir/diff_splicing_outs/DRIMSeq.out --tx2gene /home/indices/R/$indexAppendix/tx2gene.RData" && $drimseq_script --counts $outDir/KALLISTO/alignment --pdata $pData --outfile $outDir/diff_splicing_outs/DRIMSeq.out\
	 --tx2gene /home/indices/R/$indexAppendix/tx2gene.RData )

echo "[INFO] [DRIMSeq] ["`date "+%Y/%m/%d-%H:%M:%S"`"] Splicing analysis finished"$'\n'
