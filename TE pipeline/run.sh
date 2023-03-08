#!/bin/bash

SECONDS=0

#run EDTA
#mkdir output
#cd output/
#mkdir EDTA/
#cd EDTA/

#conda activate EDTA
#perl ../../../../software/EDTA/EDTA.pl --genome ../../genome/*.fa --anno 1 --threads 45
#conda deactivate

#cd ../../

cd genome/

#get fai file
samtools faidx *.fa

cd ../

cd output/
mkdir filter python
cd filter/

#only chromosomes, no scaffolds
awk '$1 == ($1+0)' ../EDTA/*.TEanno.gff3 > chromosomes.TEanno.gff3
cp chromosomes.TEanno.gff3 seqextraction.txt

cp seqextraction.txt ../python

conda activate base

#filter sequences, remove sequences < 100 bp, if same start or end only keep the longest sequence, remove helitrons
python3 ../../scripts/python/edtasectionsandfilter.py

rm seqextraction.txt 

cd ../python/

cp ../filter/edta_seqe_filter_noh.txt edta_seqe_filter_noh.txt

#obtain sequences
bedtools slop -l 20 -r 20 -i ../filter/edta_seqe_filter_noh.bed.txt -g ../../genome/*.fa.fai > edta_seqe_filter_noh20.bed.txt
bedtools getfasta -fi ../../genome/*.toplevel.fa -bed edta_seqe_filter_noh20.bed.txt -fo outputnohseq.fa.out

#run through python algorithm for tsd/tir features and locations
one=1
touch tempvalidtsdtir.txt
touch tempvalidtsdtiredta.txt
> tempvalidtsdtir.txt
> tempvalidtsdtiredta.txt
while [ "$one" -le 10 ]
do
 python3 ../../scripts/python/tsdtirfeatures.py
 cat validtsdtir.txt >> tempvalidtsdtir.txt
 cat validtsdtiredta.txt >> tempvalidtsdtiredta.txt 
 mkdir $one
 mv tsdtir.txt tsdtirseq.bed.txt validtsdtir.txt validtsdtiredta.txt $one/
 one=$(( one+1 ))
 rm -rf $one/
done

mkdir all/
cd all/

cp ../tempvalidtsdtir.txt tempvalidtsdtir.txt
cp ../tempvalidtsdtiredta.txt tempvalidtsdtiredta.txt

rm ../temp*.txt

#remove repeats
awk '{print $2}' tempvalidtsdtir.txt | sort | uniq -D | awk '!seen[$1]++' > uniqcode.txt
grep -Fwf uniqcode.txt tempvalidtsdtir.txt | sort -k4n -k5n -k8 | awk '!seen[$2]++' > rptvalidtsdtir.txt
grep -Fwf uniqcode.txt tempvalidtsdtiredta.txt | sort -k4n -k5n -k8 | awk '!seen[$2]++' > rptvalidtsdtiredta.txt

cd ../

cp all/rptvalidtsdtir.txt validtsdtir.txt
cp all/rptvalidtsdtiredta.txt validtsdtiredta.txt

rm -rf all/

conda deactivate

#match location
awk '{print $3":"$6"-"$7}' validtsdtiredta.txt > locvalidtsdtiredta.txt
awk '{print $4":"$5"-"$6}' validtsdtir.txt > locvalidtsdtir.txt
grep -Fwf locvalidtsdtir.txt locvalidtsdtiredta.txt > matchvaledta.txt
paste validtsdtir.txt locvalidtsdtir.txt > newvalidtsdtir.txt
paste validtsdtiredta.txt locvalidtsdtiredta.txt > newvalidtsdtiredta.txt
grep -Fwf matchvaledta.txt newvalidtsdtir.txt | sort -k15h -k8 | awk '!seen[$15]++' | sort -k4n -k5n > newmatchedta.txt
grep -Fwf matchvaledta.txt newvalidtsdtiredta.txt > newmatchedta2.txt
awk -F "\t" 'OFS="\t" {print $4, $5, $6}' newmatchedta.txt | sed '1d' > match.bed.txt
paste <(awk -F "\t" 'OFS="\t" {print $2, $4, $5, $6 }' newmatchedta.txt ) <(awk -F "\t" 'OFS="\t" {print $9 }' newmatchedta2.txt ) <(awk -F "\t" 'OFS="\t" {print $7, $8 }' newmatchedta.txt ) <(awk -F "\t" 'OFS="\t" {print $11, $13, $18, $19 }' newmatchedta2.txt ) <(awk -F "\t" 'OFS="\t" {print $9, $11 }' newmatchedta.txt | tr [:lower:] [:upper:] ) <(awk -F "\t" 'OFS="\t" {print $13 }' newmatchedta.txt ) > lmatchedta.txt

#rm locvalidtsdtiredta.txt locvalidtsdtir.txt matchvaledta.txt newvalidtsdtir.txt newvalidtsdtiredta.txt newmatchedta.txt newmatchedta2.txt

#match classification
awk -F "\t" 'OFS="\t" {print $1 }' lmatchedta.txt > codematch.txt
grep -Fwf codematch.txt validtsdtiredta.txt | sort -k3n -k6n | awk -F "\t" 'OFS="\t" {print $2, $14, $15 }' > classedta.txt
awk -F "\t" 'OFS="\t" {print $1, $2, $3, $4, $6, $7 }' lmatchedta.txt > classalg.txt
paste classalg.txt classedta.txt | grep 'LTR' > ltr.txt
paste classalg.txt classedta.txt | grep -v 'LTR' | grep -v 'OTH' > dna.txt
paste classalg.txt classedta.txt | grep -v 'LTR' | grep 'OTH' > oth.txt
#check DNA classification
awk '{if(substr($6,0,3)==$9) print $0}' dna.txt > rightclassdna.txt
awk -F "\t" 'OFS="\t" {if(substr($6,0,3)!=$9) print $1, $2, $3, $4, $5, "unk", $7, $8, $9}' dna.txt > corrclassdna.txt
#check LTR classification
awk '{if($5==$8) print $0}' ltr.txt > rightclassltr.txt
awk -F "\t" 'OFS="\t" {if($5!=$8) print $1, $2, $3, $4, "UNK", "unk", $7, $8, $9}' ltr.txt > corrclassltr.txt
#check OTH classification
awk -F "\t" 'OFS="\t" {print $1, $2, $3, $4, "UNK", "unk", $7, $8, $9}' oth.txt > corrclassoth.txt

cat rightclassdna.txt corrclassdna.txt rightclassltr.txt corrclassltr.txt corrclassoth.txt | sort -k2n -k3n > modclass.txt

paste <(awk -F "\t" 'OFS="\t" {print $1, $2, $3, $4, $5 }' lmatchedta.txt ) <(awk -F "\t" 'OFS="\t" {print $5, $6 }' modclass.txt ) <(awk -F "\t" 'OFS="\t" {print $8, $9, $10, $11, $12, $13, $14 }' lmatchedta.txt ) > matchedta.txt

#rm codematch.txt classedta.txt classalg.txt ltr.txt dna.txt rightclassdna.txt corrclassdna.txt rightclassltr.txt corrclassltr.txt modclass.txt oth.txt corrclassoth.txt

cd ../
mkdir blat
cd blat/

#fasta of structure matches for blat
bedtools slop -l 1 -r 0 -i ../python/match.bed.txt -g ../../genome/*.fa.fai > l1match.bed.txt
bedtools getfasta -fi ../../genome/*.toplevel.fa -bed l1match.bed.txt -fo outputl1match.fa.out
grep '^>' outputl1match.fa.out | awk 'sub(/^>/, "")' > referencel1match.txt
awk '{ print $8 }' ../python/matchedta.txt | sed '1d' > match_id.txt
paste referencel1match.txt match_id.txt > reference_match.txt
awk 'NR==FNR{a[$1]=$2;next} NF==2{$2=a[$2]; print ">" $2;next}1' FS='\t' reference_match.txt FS='>' outputl1match.fa.out > idmatch.fa.out

#pblat
pblat ../../genome/*.toplevel.fa idmatch.fa.out -threads=8 -tileSize=8 -stepSize=2 -minIdentity=95 -repMatch=1000000 blatmatch.psl

cp ../python/matchedta.txt matchedta.txt

#remove headline
sed '1,5d' blatmatch.psl > nohblatmatch.txt
#get columns of interest zero-based
awk -F "\t" 'OFS="\t" { print $10, $1, $12, $13, $11, $14, $9, $16, $17, $17-$16, ($17-$16)*100/$11, $14":"$16"-"$17 }' nohblatmatch.txt > blatmatch.txt
#remove scaffolds
awk '$6 == ($6+0)' blatmatch.txt > blatmatchc.txt

#sort
sort -k6n -k8n -k10n -k11hr blatmatchc.txt > sortblatmatchc.txt
grep -vf referencel1match.txt sortblatmatchc.txt > noqblat.txt
#grep -Fw -f referencel1match.txt sortblatmatchc.txt > queryblat.txt
awk '!seen[$12]++' noqblat.txt > 12noqblat.txt
sort -k6n -k8n -k10nr 12noqblat.txt > sort12noqblat.txt
awk '!seen[$8]++' sort12noqblat.txt > 8noqblat.txt
sort -k6n -k9n -k10nr 8noqblat.txt > sort8noqblat.txt
awk '!seen[$9]++' sort8noqblat.txt > 9noqblat.txt

#cat 1queryblat.txt 9noqblat.txt > blatfilterz.txt
awk -F "\t" 'OFS="\t" { print $1, $2, $3, $4, $5, $6, $7, $8+1, $9 }' 9noqblat.txt > blatfilter.txt
#get bed file
awk -F "\t" 'OFS="\t" { print $6, $8, $9 }' blatfilter.txt > blatmatchc.bed.txt

#add 20bp on each side of sequence, for query and blat results
bedtools slop -l 21 -r 20 -i blatmatchc.bed.txt -g ../../genome/*.fa.fai > 20blatmatchc.bed.txt
bedtools slop -l 20 -r 20 -i l1match.bed.txt -g ../../genome/*.fa.fai > locationqueryblat.txt
#get fasta of blat results (l1)
bedtools getfasta -fi ../../genome/*.toplevel.fa -bed 20blatmatchc.bed.txt -fo l1blatmatchc.fa.out
#get fasta titles of blat results
grep '^>' l1blatmatchc.fa.out | awk 'sub(/^>/, "")' > referencel1mc.txt
#fasta like name of query sequences and blat results, one-based
bedtools slop -l 20 -r 20 -i ../python/match.bed.txt -g ../../genome/*.fa.fai > 20locationqueryblat.txt
awk '{ $1=$1":"$2"-"$3;$2="";$3="" } 1' 20locationqueryblat.txt > locquery.txt
bedtools slop -l 20 -r 20 -i blatmatchc.bed.txt -g ../../genome/*.fa.fai > 20oblatmatchc.bed.txt
awk '{ $1=$1":"$2"-"$3;$2="";$3="" } 1' 20oblatmatchc.bed.txt > locblatmc.txt
#zero-based fasta names next two one-based fasta names
paste referencel1mc.txt locblatmc.txt > referencematch.txt
awk 'NR==FNR{a[$1]=$2;next} NF==2{$2=a[$2]; print ">" $2;next}1' FS='\t' referencematch.txt FS='>' l1blatmatchc.fa.out > blatmatch.fa.out
#join fasta names and blat results
paste locblatmc.txt blatfilter.txt > blatmatchf.txt

conda activate base

#find blat tsd and tir features
python3 ../../scripts/python/blatfeatures.py

conda deactivate

awk -F "\t" 'OFS="\t" {print $2, "algorithm", $6"/"$7, $3, $4, ".", $5, ".", "TSD="$10";TIR="$11";length="$4-$3+1 }' blatfeatures.txt | sed '1d' > blatfeaturesnoq.txt
awk -F "\t" 'OFS="\t" {print $2, "algorithm", $6"/"$7, $3, $4, ".", $5, ".", "TSD="$12";TIR="$13";length="$4-$3+1 }' matchedta.txt > blatfeaturesq.txt

cat blatfeaturesnoq.txt blatfeaturesq.txt > blatfeaturesall.txt

#rm l1match.bed.txt outputl1match.fa.out referencel1match.txt match_id.txt reference_match.txt

#rm locblatmc.txt blatfilter.txt referencematch.txt locblatmc.txt 20oblatmatchc.bed.txt locquery.txt 20locationqueryblat.txt referencel1mc.txt l1blatmatchc.fa.out locationqueryblat.txt 20blatmatchc.bed.txt blatmatchc.bed.txt blatfilter.txt matnoqblat.txt refbestm.txt bestmnoqblat.txt noqblat.txt sortblatmatchc.txt perblatmatchc.txt 100blatmatchc.txt blatmatchc.txt blatmatch.txt nohblatmatch.txt

cd ../../
mkdir results
cd results/

sort -k1n -k4n -k5n ../output/blat/blatfeaturesall.txt > TEfinal.gff

awk -F "\t" 'OFS="\t" { print $1, $4-1, $5 }' TEfinal.gff > TEfinal.bed
awk -F "\t" 'OFS="\t" { print $1":"$4"-"$5"_"$3"_"$7"_"$9 }' TEfinal.gff > TEfinalref.txt
bedtools getfasta -fi ../genome/*.toplevel.fa -bed TEfinal.bed -fo TEfinall1.fa
grep '^>' TEfinall1.fa | awk 'sub(/^>/, "")' > referenceTEfinalfa.txt
paste referenceTEfinalfa.txt TEfinalref.txt > referencefa.txt
awk 'NR==FNR{a[$1]=$2;next} NF==2{$2=a[$2]; print ">" $2;next}1' FS='\t' referencefa.txt FS='>' TEfinall1.fa > TEfinal.fa

rm referencefa.txt referenceTEfinalfa.txt TEfinalref.txt TEfinall1.fa

sec=$SECONDS
printf 'Run time: %dh:%dm:%ds\n' $((sec/3600)) $((sec%3600/60)) $((sec%60))

