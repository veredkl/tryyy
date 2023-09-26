$HOSTNAME = ""
params.outdir = 'results'  


if (!params.inputparam1){params.inputparam1 = ""} 
if (!params.inputparam){params.inputparam = ""} 

Channel.fromPath(params.inputparam1, type: 'any').map{ file -> tuple(file.baseName, file) }.set{g_1_0_g0_21}
Channel.value(params.inputparam).into{g_2_0_g0_9;g_2_1_g0_23;g_2_0_g0_25}

//* params.run_Remove_Multimappers_with_Samtools =  "no"  //* @dropdown @options:"yes","no" @show_settings:"Remove_Multimappers"

if (!(params.run_Remove_Multimappers_with_Samtools == "yes")){
g_1_0_g0_21.set{g0_21_mapped_reads00_g0_22}
} else {


process ChIP_Module_Remove_Multimappers_with_Samtools {

input:
 set val(name), file(bam) from g_1_0_g0_21

output:
 set val(name), file("bam/${name}.bam")  into g0_21_mapped_reads00_g0_22

when:
params.run_Remove_Multimappers_with_Samtools == "yes" 

script:
MAPQ_quality = params.ChIP_Module_Remove_Multimappers_with_Samtools.MAPQ_quality
"""
mkdir bam
samtools view -hb -q ${MAPQ_quality} ${bam} > ${name}_unique.bam
mv ${name}_unique.bam bam/${name}.bam
"""

}
}



//* autofill
if ($HOSTNAME == "default"){
    $CPU  = 1
    $MEMORY = 32
}
//* platform
if ($HOSTNAME == "ghpcc06.umassrc.org"){
    $TIME = 240
    $CPU  = 1
    $MEMORY = 32
    $QUEUE = "short"
}
//* platform
//* autofill
if (!((params.run_Remove_Multimappers_with_Picard && (params.run_Remove_Multimappers_with_Picard == "yes")) || !params.run_Remove_Multimappers_with_Picard)){
g0_21_mapped_reads00_g0_22.set{g0_22_mapped_reads01_g0_25}
g0_22_publish11 = Channel.empty()
g0_22_log_file20_g0_23 = Channel.empty()
} else {


process ChIP_Module_Picard_MarkDuplicates {

input:
 set val(name), file(bam) from g0_21_mapped_reads00_g0_22

output:
 set val(name), file("bam/${name}.bam")  into g0_22_mapped_reads01_g0_25
 set val(name), file("${name}*")  into g0_22_publish11
 file "*_duplicates_stats.log"  into g0_22_log_file20_g0_23

when:
(params.run_Remove_Multimappers_with_Picard && (params.run_Remove_Multimappers_with_Picard == "yes")) || !params.run_Remove_Multimappers_with_Picard     

script:
"""
mkdir bam
picard MarkDuplicates OUTPUT=${name}_dedup.bam METRICS_FILE=${name}_picard_PCR_duplicates.log  VALIDATION_STRINGENCY=LENIENT REMOVE_DUPLICATES=false INPUT=${bam} > ${name}_picard.log 

#get duplicates stats (read the sam flags)
samtools flagstat ${name}_dedup.bam > ${name}@Reads@${name}_duplicates_stats.log
#remove alignments marked as duplicates
samtools view -b -F 0x400 ${name}_dedup.bam > ${name}_dedup.bam.x_dup
#sort deduplicated files by chrom pos
samtools sort -o ${name}_sorted.dedup.bam ${name}_dedup.bam.x_dup 
mv ${name}_sorted.dedup.bam bam/${name}.bam
#get properly paired log stats after dedup
echo "##After Deduplication##" >> ${name}@Reads@${name}_duplicates_stats.log
samtools flagstat bam/${name}.bam >> ${name}@Reads@${name}_duplicates_stats.log
"""
}
}


macs2_callpeak_parameters = params.ChIP_Module_ChIP_Prep.macs2_callpeak_parameters
peak_calling_type = params.ChIP_Module_ChIP_Prep.peak_calling_type
band_width = params.ChIP_Module_ChIP_Prep.band_width
bedtoolsCoverage_Parameters = params.ChIP_Module_ChIP_Prep.bedtoolsCoverage_Parameters
compare_Custom_Bed = params.ChIP_Module_ChIP_Prep.compare_Custom_Bed
output_prefix = params.ChIP_Module_ChIP_Prep.output_prefix
sample_prefix = params.ChIP_Module_ChIP_Prep.sample_prefix
input_prefix = params.ChIP_Module_ChIP_Prep.input_prefix
//* @array:{output_prefix,sample_prefix,input_prefix} @multicolumn:{output_prefix,sample_prefix,input_prefix},{macs2_callpeak_parameters,peak_calling_type,band_width,bedtoolsCoverage_Parameters}
samplehash = [:]
inputhash = [:]
output_prefix.eachWithIndex { key, i -> inputhash[key] = input_prefix[i] }
output_prefix.eachWithIndex { key, i -> samplehash[key] = sample_prefix[i] }

process ChIP_Module_ChIP_Prep {

input:
 val mate from g_2_0_g0_25
 set val(name), file(bam) from g0_22_mapped_reads01_g0_25

output:
 file "bam/*.bam"  into g0_25_bam_file01_g0_9
 val output_prefix  into g0_25_name12_g0_9

when:
(params.run_ChIP_MACS2 && (params.run_ChIP_MACS2 == "yes")) || !params.run_ChIP_MACS2

script:
"""
mkdir -p bam
mv ${bam} bam/${name}.bam
"""
}


process ChIP_Module_ChIP_MACS2 {

input:
 val mate from g_2_0_g0_9
 file bam from g0_25_bam_file01_g0_9.collect()
 val name from g0_25_name12_g0_9.unique().flatten()

output:
 val compare_bed  into g0_9_compare_bed00_g0_27
 file "*${peak_calling_type}Peak"  into g0_9_bed10_g0_10
 set val(name), file("bam/*.bam")  into g0_9_bam_file21_g0_10, g0_9_bam_file22_g0_27
 file "${name}*"  into g0_9_resultsdir33
 val name  into g0_9_name44

script:
genomeSizeText = ""
if (params.genome_build.contains("mouse")){
    genomeSizeText = "-g mm"
} else if (params.genome_build.contains("human")){
    genomeSizeText = "-g hs"
}

if (peak_calling_type == "narrow"){
    peakcallingType = ""
} else if (peak_calling_type == "broad"){
    peakcallingType = "--broad"
}

compare_bed = "merged.bed"
compare_Custom_Bed = compare_Custom_Bed.trim();
if (compare_Custom_Bed != ""){
    compare_bed = compare_Custom_Bed
}
inputsList = inputhash[name] 
samplesList = samplehash[name]

"""
echo ${samplesList}
echo ${inputsList}
echo $name
mkdir -p bam

#samplesList
samplesList="\$(echo -e "${samplesList}" | tr -d '[:space:]')" 
IFS=',' read -ra eachSampleAr <<< "\${samplesList}"
numSamples=\${#eachSampleAr[@]}
eachSampleArBam=( "\${eachSampleAr[@]/%/.bam }" )
sample_set=\${eachSampleArBam[@]}
bam_set=\${eachSampleArBam[@]}

#inputsList
input_set=""
inputsList="\$(echo -e "${inputsList}" | tr -d '[:space:]')" 
if [ "\${inputsList}" != "" ]; then
    IFS=',' read -ra eachInputAr <<< "\${inputsList}"
    eachInputArbam=( "\${eachInputAr[@]/%/.bam }" )
    input_set="-c \${eachInputArbam[@]}" 
    
fi
echo \${eachSampleArBam[@]}

macs2 callpeak --bw ${band_width} -t \${sample_set} \${input_set} -n ${name} ${genomeSizeText} ${macs2_callpeak_parameters} ${peakcallingType}

#bam files
if [ "\$numSamples" -gt "1" ]; then
    samtools merge bam/${name}.bam \$bam_set
else 
    rsync -a  \$bam_set bam/${name}.bam
fi

"""
}

//* params.run_Scripture =  "no"  //* @dropdown @options:"yes","no" @show_settings:"Scripture_peakrescore"
//* params.peakrescore_path =  ""  //* @input
//* params.peakrescore_class_path =  ""  //* @input

if (!(params.run_Scripture == "yes")){
g0_9_bed10_g0_10.set{g0_10_bed00_g0_26}
} else {


process ChIP_Module_Scripture_peakrescore {

input:
 file bed from g0_9_bed10_g0_10
 set val(name), file(bam) from g0_9_bam_file21_g0_10

output:
 file "${name}_trim.bed"  into g0_10_bed00_g0_26

when:
params.run_Scripture == "yes"

script:
window = params.ChIP_Module_Scripture_peakrescore.window
trimFraction = params.ChIP_Module_Scripture_peakrescore.trimFraction
windowText = (window.toString() != "") ? "-window ${window}" : ""
trimFractionText = (trimFraction.toString() != "") ? "-trimFraction ${trimFraction}" : ""
"""
samtools index ${bam}
cat ${bed} | awk '{print \$1"\t"\$2"\t"\$3"\t"\$4"\t"\$5}' > ${name}_clean 
java -cp ${params.peakrescore_path}:${params.peakrescore_class_path} peaks.PeakTrim -task trimByFractionOfScale -in ${name}_clean -libAlignment ${bam}  $windowText $trimFractionText -out ${name}_trim.bed 
"""

}
}



process ChIP_Module_bed_merge {

input:
 file bed from g0_10_bed00_g0_26.collect()

output:
 file "merged.bed"  into g0_26_bed01_g0_27

"""
 cat ${bed} | cut -f -6 | bedtools sort -i stdin | bedtools slop -i stdin -b 100 -g ${params.genome_sizes} | bedtools merge -i stdin | awk '{print \$0"\t"\$1"_"\$2"_"\$3}' > merged.bed

"""
}



process ChIP_Module_bedtools_coverage {

input:
 val compare_bed from g0_9_compare_bed00_g0_27
 file bed from g0_26_bed01_g0_27
 set val(name), file(bam) from g0_9_bam_file22_g0_27

output:
 file "*.sum.txt"  into g0_27_outputFileTxt00_g0_13

script:
bedtoolsCoverage_Parameters = params.ChIP_Module_bedtools_coverage.bedtoolsCoverage_Parameters
bedtoolsIntersect_Parameters = params.ChIP_Module_bedtools_coverage.bedtoolsIntersect_Parameters
"""
echo ${compare_bed}
if [ -s "${compare_bed}" ]; then 
    echo " bed file exists and is not empty "
        samtools view -H ${name}.bam | grep -P "@SQ\\tSN:" | sed 's/@SQ\\tSN://' | sed 's/\\tLN:/\\t/' > ${name}_chroms
        samtools sort -T ${name} -o ${name}_sorted.bam ${name}.bam
        bedtools intersect -abam ${name}_sorted.bam -b ${compare_bed} > temp_${name}.bam
        bedtools sort -faidx ${name}_chroms -i ${compare_bed}  | bedtools coverage ${bedtoolsCoverage_Parameters} -a stdin -b temp_${name}.bam  > temp_${name}.bed
        # 'The number of features in B that overlapped the A interval' multiplied by 'fraction of bases in A that had non-zero coverage from features in B'.
        awk '{\$NF=\$(NF-3)*\$NF;print }' OFS="\\t" temp_${name}.bed | grep -v all > temp_${name}_hist.bed
        l=`awk '{print NF}' temp_${name}_hist.bed | head -1 | awk '{print \$1-4}'`
        k=`awk '{print NF}' temp_${name}_hist.bed | head -1`
        bedtools groupby -i temp_${name}_hist.bed -g 1-\$l -c \$k -o sum > ${name}.sum.txt
        #rm -rf temp_*

else
  echo " bed file does not exist, or is empty "
  touch ${name}_empty.sum.txt
fi
"""

}


process ChIP_Module_ATAC_CHIP_summary {

input:
 file file from g0_27_outputFileTxt00_g0_13.collect()

output:
 file "*.tsv"  into g0_13_outputFile00

shell:
'''
#!/usr/bin/env perl

my $indir = $ENV{'PWD'};

opendir D, $indir or die "Could not open $indir\n";
my @alndirs = sort { $a cmp $b } grep /.txt/, readdir(D);
closedir D;
    
my @a=();
my %b=();
my %c=();
my $i=0;
foreach my $d (@alndirs){ 
    my $file = "${indir}/$d";
    print $d."\n";
    my $libname=$d;
    $libname=~s/\\.sum\\.txt//;
    print $libname."\n";
    $i++;
    $a[$i]=$libname;
    open IN,"${indir}/$d";
    $_=<IN>;
    while(<IN>)
    {
        my @v=split; 
        $b{$v[3]}{$i}=$v[4];
    }
    close IN;
}
my $outfile="${indir}/"."sum_counts.tsv";
open OUT, ">$outfile";
print OUT "Feature";

for(my $j=1;$j<=$i;$j++) {
    print OUT "\t$a[$j]";
}
print OUT "\n";
    
foreach my $key (keys %b){
    print OUT "$key";
    for(my $j=1;$j<=$i;$j++){
        print OUT "\t$b{$key}{$j}";
    }
    print OUT "\n";
}
close OUT;
'''
}


process ChIP_Module_Picard_Deduplication_Summary {

input:
 file flagstat from g0_22_log_file20_g0_23.collect()
 val mate from g_2_1_g0_23

output:
 file "deduplication_summary.tsv"  into g0_23_outputFileTSV00

errorStrategy 'retry'
maxRetries 2

shell:
'''
#!/usr/bin/env perl
use List::Util qw[min max];
use strict;
use File::Basename;
use Getopt::Long;
use Pod::Usage;
use Data::Dumper;

my @header;
my %all_files;
my %tsv;
my %headerHash;
my %headerText;

my $i=0;
chomp(my $contents = `ls *_duplicates_stats.log`);
my @files = split(/[\\n]+/, $contents);
foreach my $file (@files){
    $i++;
    $file=~/(.*)@(.*)@(.*)_duplicates_stats\\.log/;
    my $mapOrder = int($1); 
    my $mapper = $2; #mapped element 
    my $name = $3; ##sample name
    push(@header, $mapper) unless grep{$_ eq $mapper} @header; 
        
    # my $duplicates;
    my $aligned;
    my $dedup; #aligned reads after dedup
    my $percent=0;
    if ("!{mate}" eq "pair" ){
        #first flagstat belongs to first bam file
        chomp($aligned = `cat $file | grep 'properly paired (' | sed -n 1p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
        #second flagstat belongs to dedup bam file
        chomp($dedup = `cat $file | grep 'properly paired (' | sed -n 2p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    } else {
        chomp($aligned = `cat $file | grep 'mapped (' | sed -n 1p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
        chomp($dedup = `cat $file | grep 'mapped (' | sed -n 2p | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    }
    # chomp($duplicates = `cat $file | grep 'duplicates' | awk '{sum+=\\$1+\\$3} END {print sum}'`);
    # $dedup = int($aligned) - int($duplicates);
    if ("!{mate}" eq "pair" ){
       $dedup = int($dedup/2);
       $aligned = int($aligned/2);
    } 
    $percent = "0.00";
    if (int($aligned)  > 0 ){
       $percent = sprintf("%.2f", ($aligned-$dedup)/$aligned*100); 
    } 
    $tsv{$name}{$mapper}=[$aligned,$dedup,"$percent%"];
    $headerHash{$mapOrder}=$mapper;
    $headerText{$mapOrder}=["$mapper (Before Dedup)", "$mapper (After Dedup)", "$mapper (Duplication Ratio %)"];
}

my @mapOrderArray = ( keys %headerHash );
my @sortedOrderArray = sort { $a <=> $b } @mapOrderArray;

my $summary = "deduplication_summary.tsv";
open(OUT, ">$summary");
print OUT "Sample\\t";
my @headArr = ();
for my $mapOrder (@sortedOrderArray) {
    push (@headArr, @{$headerText{$mapOrder}});
}
my $headArrAll = join("\\t", @headArr);
print OUT "$headArrAll\\n";

foreach my $name (keys %tsv){
    my @rowArr = ();
    for my $mapOrder (@sortedOrderArray) {
        push (@rowArr, @{$tsv{$name}{$headerHash{$mapOrder}}});
    }
    my $rowArrAll = join("\\t", @rowArr);
    print OUT "$name\\t$rowArrAll\\n";
}
close(OUT);
'''
}


workflow.onComplete {
println "##Pipeline execution summary##"
println "---------------------------"
println "##Completed at: $workflow.complete"
println "##Duration: ${workflow.duration}"
println "##Success: ${workflow.success ? 'OK' : 'failed' }"
println "##Exit status: ${workflow.exitStatus}"
}
