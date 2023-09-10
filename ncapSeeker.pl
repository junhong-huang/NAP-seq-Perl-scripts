#!/usr/bin/perl -w

BEGIN {
    unshift( @INC, "/public/home/huangjh/perlBioTools/packages" );
}
use Getopt::Long;
use strict;
use commonTools;
use bioUtils;
use fastTools;
use alignTools;
use samTools;
use annoTools;
use evolutionTools;
use Math::Complex;
use XML::Simple;
use File::Basename;

$| = 1;

my $configureXmlFile = "xml";
my $read1File        = "";
my $read2File        = "";
my $outputDir        = "./output";
my $prefix           = "prefix";
my $aligner          = "bowtie2";
my $platform         = "pairedEnd";
my $readLen          = 150;
my $seqFormat        = "fastq";
my $help             = 0;
my $org              = "hg38";
my $extLen = 20;

if ( scalar(@ARGV) < 4 ) {
    usage();
}

my $optLong = GetOptions(
    "conf=s"     => \$configureXmlFile,    # configure file
    "read1=s"    => \$read1File,           # read1File file
    "read2=s"    => \$read2File,           # read2File file
    "output=s"   => \$outputDir,           # outputDir
    "prefix=s"   => \$prefix,              # prefix
    "aligner=s"  => \$aligner,             # aligner
    "platform=s" => \$platform,            # platform
    "readLen=s"  => \$readLen,             # readLen
    "format=s"   => \$seqFormat,           # sequence format
    "org=s"      => \$org,                 # org
    "help|h"     => \$help
) or usage();

&usage() if $help;

### configure files
my $xmlConfHash = XMLin($configureXmlFile);
my $confHash = XMLin($configureXmlFile);

sub usage
### usage function
{
    print
        "Usage: perl $0 [options] --conf <configure file> --aligner <bowtie2> --output <output file> --prefix <file prefix> 
        --platform <singleEnd> --readLen <150> --read1 <read1 file> --read2 <read2 file> --org <organism, e.g. hg38> \n";
    print "--conf      conf file\n";
    print "--read1     read1  file\n";
    print "--read2     read2  file\n";
    print "--output    output dir name[default=output]\n";
    print "--prefix    prefix for file\n";
    print
        "--aligner     aligner tool[default=bowtie2], bowtie or bowtie2 or star\n";
    print
        "--platform    platform [default=singleEnd], pairedEnd or singleEnd\n";
    print "--readLen   read length [default=101]\n";
    print "--format    seq format [default=fastq], fastq or fasta\n";
    print "--org       organism [default=hg38]\n";
    print "--help      help information\n";
    exit;
}

if ( !( -e $outputDir ) ) {
    warn "ceate output dir\n";
    &executeCommand("mkdir $outputDir");
}

our $demo=0;

my $SHAPE=" 
cd /public/home/huangjh/NAPseq/SHAPE-MaP && \
nohup perl /public/home/huangjh/perlBioTools/ncapSeeker/doNcapSeeker.pl \
--samp merge_HepG2_NAPseq_SHAPE_samples.txt &
";

#### TGS
my $NcapNanopore=" 
cd /public/home/huangjh/NAPseq/nanopore_3NGS && \
nohup perl /public/home/huangjh/perlBioTools/ncapSeeker/doNcapSeeker.pl \
--samp human_NcapNanopore_HepG2_samples.txt &
";

my $NcapNanopore=" 
cd /public/home/huangjh/NAPseq/nanopore_3NGS && \
nohup perl /public/home/huangjh/perlBioTools/ncapSeeker/doNcapSeeker.pl \
--samp mouse_NcapNanopore_C2C12_samples.txt &
";

my $HepG2_rep=" 
cd /public/home/huangjh/NAPseq/human && \
nohup perl /public/home/huangjh/perlBioTools/ncapSeeker/doNcapSeeker.pl \
--samp HepG2_rep1_samples.txt &
";

my $C2C12=" 
cd /public/home/huangjh/NAPseq/mouse && \
nohup perl /public/home/huangjh/perlBioTools/ncapSeeker/doNcapSeeker.pl \
--samp C2C12_all_samples.txt &
";



$prefix.="_splice" if ($outputDir~~/nanopore/i);
# global parameters
my $species="human";
$species="mouse" if ($org~~/mm10/);
our $minPair  = 10;
our $maxMfe   = -10.0;
our $minScore = 5;

#1.计算多少read含有特定接头,输出统计信息(文件.startEndTagSummary.txt)
#my $readSummaryFile = &computeStartEndTag($xmlConfHash, $read1File, $outputDir, $prefix);

#2.比对
my $bamFile = $outputDir . "/" . $prefix . ".STAR.Aligned.out.bam";
$bamFile = $outputDir . "/" . $prefix . ".minimap2.sorted.bam" if ($outputDir~~/nanopore/i);
# if ( (-s $bamFile)/1024 < 100 ){
#   my ($bamFile, $unmappedRead1File, $unmappedRead2File) = &alignReadToGenome(
#     $read1File, $read2File, $xmlConfHash, $outputDir,
#     $prefix,    $aligner,   $platform,    $seqFormat
#     );  
# }

# my $intronScanFile = &runIntronScan( $xmlConfHash, $bamFile, $outputDir, $prefix );
# print "$bamFile < 100K\n";

# my $bedAnnotatorSummaryFile = &annotateResultsSummary($xmlConfHash, $bamFile);
my $bedAnnotatorSummaryFile=$outputDir . "/" . $prefix . ".STAR.Aligned.out.bam.bedAnnotator.summary.txt";
# &pieplot_common( $xmlConfHash, $bedAnnotatorSummaryFile, $outputDir, 0, 100);


# ##### intronScan #######
# ##### intronScan #######
# ##### intronScan #######
# my $intronScanFile = &runIntronScan( $xmlConfHash, $bamFile, $outputDir, $prefix );
# my $intronScanFile = $outputDir . "/" . $prefix . ".intronScan";
# my $mergeFile = &mergeIntron( $xmlConfHash, $intronScanFile, $outputDir, $prefix );


##### rriScan #######
##### rriScan #######
##### rriScan #######
# my $bamFile = $outputDir . "/" . $prefix . ".STAR.Aligned.sortedByCoord.out.bam";
# my $junFile = $outputDir . "/" . $prefix . ".STAR.Chimeric.out.junction";
# my $rriScanResults = $outputDir . "/" . $prefix . ".rriScan";
# if ( (-s $bamFile)/1024 < 1000 ){
#     my ($bamFile, $unmappedRead1File, $unmappedRead2File) = &alignReadToGenome_rriScan(
#         $read1File, $read2File, $xmlConfHash, $outputDir,
#         $prefix,    $aligner,   $platform,    $seqFormat
#     );
#     $rriScanResults = &runRriScan( $xmlConfHash, $bamFile, $outputDir, $prefix, $junFile );
# }
# my $annoteResults = &annotateChimeras( $xmlConfHash, $rriScanResults, $outputDir, $prefix );
# my $filterNum     = &filterCandidates( $xmlConfHash, $annoteResults, $outputDir, $prefix );
# my $ratioFile     = &getChimericRatio( $xmlConfHash, $filterNum, $outputDir, $prefix );



#3.调用NcapSeeker（什么来的），输入$bamFile->输出$ncapFile(文件.ncapSeeker)；
# my $ncapFile = &runNcapSeeker( $xmlConfHash, $bamFile, $outputDir, $prefix );


# my $ncapFile= $outputDir . "/" . $prefix . ".ncapSeeker";
# #4.调用bedConservation（什么来的），输入$ncapFile；->输出$evoConValFile（文件.ncapSeeker.bedConservation;增加一列conservation值）
# my $evoConValFile = &getResultConVals( $xmlConfHash, $ncapFile);

# #5.调用bedAnnotator，输入$evoConValFile->输出$annoteResultFile(文件.ncapSeeker.bedConservation.bedAnnotator.txt；增加几列基因注释)
# my $annoteResultFile = &annotateResults( $xmlConfHash, $evoConValFile);


# my $sortrdBamFile = &generateBamIndex($bamFile);
# my $sortrdBamFile = $outputDir . "/" . $prefix . ".STAR.Aligned.out.sorted.bam";

# my $plot = &plotPofile( $bamFile, $xmlConfHash, "$outputDir", $prefix, $readLen, "pairedEnd", 0 );
# my $wiggle = &wigToBigwig( $outputDir, $prefix, $xmlConfHash );

# # if ($demo==1){
# #     my $bam_demo="/public/home/huangjh/NAPseq/human/mergeBamExpression/U15B.intersectBam.bam";
# #     my $plot = &plotPofile( $bam_demo, $xmlConfHash, "$outputDir/demo", $prefix, $readLen, "pairedEnd", 0, $primer );
# # }

#####6.调用ncapToWig，输入$bamFile->输出$wig
# my $plot    = &plotPofile_old( $xmlConfHash, $bamFile, $outputDir, $prefix );
# #####7.调用wigToBigwig，输入$wig->输出$bigwig，最后rm -f $wigFile；画图
# my $wiggle  = &wigToBigwig_nap( $xmlConfHash, $outputDir, $prefix );



#7.1 生成distSiteResult文件夹
my $distSiteOutDir = $outputDir . "/" . "distSiteResults";
# if ( !( -e $distSiteOutDir ) ) {
#     executeCommand("mkdir $distSiteOutDir");
# }

#8.调用C程序annotateSites，输入"$四个bw文件"和"$剪接位点文件"->输出距离文件，R画图（transcriptEndSite_sncRNA.EndToStart.annotateSites.txt）
# #剪接文件（转录起始/结束；外显子起始/结束；剪接5P/3P）
# &annotateTagSite( $xmlConfHash, $outputDir, $prefix, $distSiteOutDir );

#9.调用bedAnnotator，输入bam文件，输出read的转录组分布
#&annotateTagDistribution( $xmlConfHash, $bamFile, $outputDir, $prefix );


#### fetch exon end from single sample
# my $exonEndDir = $outputDir . "/" . "detachedExon";
# if ( !( -e $exonEndDir ) ) {
#     &executeCommand("mkdir -p $exonEndDir");
# }
# my $exonEndDistFile = $outputDir . "/distSiteResults/" . $prefix . ".$species". "_exonEndSite_protein_coding_rmSncRNA_rmOverlapEnd.EndToEnd.annotateSites.txt";
# executeCommand("perl /public/home/huangjh/perlBioTools/canSeq/fetchExonEnds.pl --org $org --end $exonEndDistFile --outDir $exonEndDir --prefix $prefix");

# our $annotateSitesToBed="$outputDir/distSiteResults/annotateSitesToBed";
# my $annotateSiteFile_exonEndToEnd = "$outputDir/distSiteResults/$prefix.human_exonEndSite_protein_coding_rmSncRNA_rmOverlapEnd.EndToEnd.annotateSites.txt";
# my ( $detectSite, $unDetectSite ) = &annotateSitesToBed($confHash, $annotateSiteFile_exonEndToEnd, "$annotateSitesToBed/1.exonEnd");
# plotMetaGena($confHash, $detectSite, "$annotateSitesToBed/1.exonEnd");


our $annotateSitesToBed="$outputDir/distSiteResults/annotateSitesToBed";
my $annotateSiteFile_exonEndToEnd = "$outputDir/distSiteResults/$prefix.human_transcriptEndSite_sncRNA.EndToEnd.annotateSites.txt";
my $minVal=0;
# my ( $detectSite, $unDetectSite ) = &annotateSitesToBed($confHash, $annotateSiteFile_exonEndToEnd, "$annotateSitesToBed/sncRNA", $minVal);

my $annotateSiteFile_exonStartToStart = "$outputDir/distSiteResults/$prefix.human_exonStartSite_protein_coding_rmSncRNA_rmOverlapEnd.StartToStart.annotateSites.txt";
# my ( $detectSite, $unDetectSite ) = &annotateSitesToBed($confHash, $annotateSiteFile_exonStartToStart, "$annotateSitesToBed/exonStartToStart", $minVal);

my $snoRNA="/public/home/huangjh/genome/human/hg38/annotations/snoRNA/hg38_snoRNABase.bed6";
my $snoRNAcap="/public/home/huangjh/genome/human/hg38/annotations/snoRNA/hg38_snoRNABase.bed6.cap";
my $snoRNA_U15B="/public/home/huangjh/NAPseq/human/mergeBamExpression/U15B.bed";

my $NGSnapRNA="/public/home/huangjh/NAPseq/napBase/final_results/merge/hg38.napRNA.bed";
# my $bamExpression = &runBamExpression( $xmlConfHash, $NGSnapRNA, $bamFile, $outputDir, "$prefix.final_napRNAs" );
# if ($demo==1){
#     $snoRNA="/public/home/huangjh/NAPseq/human/mergeBamExpression/U15B.bed";
#     $bamFile="/public/home/huangjh/NAPseq/human/mergeBamExpression/U15B.intersectBam.bam";
#     my $bamExpression = &runBamExpression( $xmlConfHash, $snoRNA, $bamFile, "$outputDir/demo", "$prefix.U15B" );
# }

if (1){
    foreach my $inBedFile (glob("/public/home/huangjh/NAPseq/mouse/diffExpress/mm10_refSeq_sorted_geneSymbols.bed")){
        my $suffix=(split/\./,basename($inBedFile))[0];
        &runBamExpression( $xmlConfHash, $inBedFile, $bamFile, $outputDir, "$prefix.refSeq" );
    }
}

#### diff_nap_correlation_with_host
if (0){
    if ($outputDir~~/RNAseq/ && $outputDir~~/hg38/){
        my $bed_region_1="/public/home/huangjh/NAPseq/human/diffExpress/diff_nap_correlation_with_host/HepG2_cell_lines_DEG/intersect_hostGene/final_cellDEG_hostGene.bed";
        my $bed_region_2="/public/home/huangjh/NAPseq/human/diffExpress/diff_nap_correlation_with_host/HepG2_stress_response_DEG/intersect_hostGene/final_stressDEG_hostGene.bed";
        my $bamExpression_1 = &runBamExpression( $xmlConfHash, $bed_region_1, $bamFile, $outputDir, "$prefix.cellDEG" );
        my $bamExpression_2 = &runBamExpression( $xmlConfHash, $bed_region_2, $bamFile, $outputDir, "$prefix.stressDEG" );   
    } 
    if ($outputDir!~/RNAseq/ && $outputDir~~/hg38/){
        my $bed_region_1="/public/home/huangjh/NAPseq/human/diffExpress/diff_nap_correlation_with_host/HepG2_cell_lines_DEG/intersect_hostGene/final_cellDEG_napRNA.bed";
        my $bed_region_2="/public/home/huangjh/NAPseq/human/diffExpress/diff_nap_correlation_with_host/HepG2_stress_response_DEG/intersect_hostGene/final_stressDEG_napRNA.bed";
        my $bamExpression_1 = &runBamExpression( $xmlConfHash, $bed_region_1, $bamFile, $outputDir, "$prefix.cellDEG" );
        my $bamExpression_2 = &runBamExpression( $xmlConfHash, $bed_region_2, $bamFile, $outputDir, "$prefix.stressDEG" );   
    } 

    if ($outputDir~~/RNAseq/ && $outputDir~~/mm10/){
        my $bed_region_1="/public/home/huangjh/NAPseq/mouse/diffExpress/diff_nap_correlation_with_host/intersect_hostGene/final_statusDEG_hostGene.bed";
        my $bamExpression_1 = &runBamExpression( $xmlConfHash, $bed_region_1, $bamFile, $outputDir, "$prefix.stageDEG" );
    } 
    if ($outputDir!~/RNAseq/ && $outputDir~~/mm10/){
        my $bed_region_1="/public/home/huangjh/NAPseq/mouse/diffExpress/diff_nap_correlation_with_host/intersect_hostGene/final_statusDEG_napRNA.bed";
        my $bamExpression_1 = &runBamExpression( $xmlConfHash, $bed_region_1, $bamFile, $outputDir, "$prefix.stageDEG" );
    }  
}


#10
# my $startEndFile = $outputDir . "/" . $prefix . ".STAR.Aligned.out.bam.startEndSeeker";
#my $startEndFile = &runStartEndSeeker( $bamFile, $platform, $readLen, $xmlConfHash );



# #11
# #将注释的位点分出来，留下未注释位点
# #下游分析文件$rmAnnoSitesResultFile
# my $rmAnnoSitesResultFile = &annotateKnownSites( $xmlConfHash, $startEndFile );
# my $rmCutSitesResultFile = &annotateCutSites( $xmlConfHash, $startEndFile );

# #有这文件可从这步开始
# my $siteFileSuffix = ".bed.bedAnnotator.FILTERknownSite";
# my $cuttingSitesResultFile = "./diffExpSites/" . "$prefix$siteFileSuffix";
# my $rmAnnoSitesResultFile = "$outputDir/$prefix.STAR.Aligned.out.bam.startEndSeeker.KnownSitesAnnotator.rmAnnoSites.txt";
# my $cuttingSitesResultFile = "$outputDir/$prefix.STAR.Aligned.out.bam.startEndSeeker.CutSitesAnnotator.rmAnnoCutSites.merge.filterCutSites.bed";

# #进一步筛选出切割位点
# my $mergeSiteFile = "$outputDir/$prefix.STAR.Aligned.out.bam.startEndSeeker.CutSitesAnnotator.rmAnnoCutSites.merge.bed";
# executeCommand("bedtools merge -d 0 -nms -s -n -i $rmCutSitesResultFile > $mergeSiteFile");
# my $cuttingSitesResultFile = &cuttingSiteFilter( $mergeSiteFile );

# # 12
# # 未注释位点进行基因组信息的注释，看看未注释位点落在哪里
# &annotateResults($xmlConfHash, $cuttingSitesResultFile);
# my $bedAnnotatorSummaryFile = &annotateResultsSummary($xmlConfHash, $cuttingSitesResultFile);
# &mySummaryToPlotPie($bedAnnotatorSummaryFile, "$prefix.cuttingSites");

#13 motif
# &motifFinder( $cuttingSitesResultFile, $outputDir, "$prefix$siteFileSuffix", $org, $extLen );

# # #14 RBP
# &executeCommand("perl /public/home/huangjh/NAPseq/codes/ncapSeeker/stableIntrons/fetchRbpIntron.pl --outDir $outputDir --org $org --stable $cuttingSitesResultFile --prefix $prefix$siteFileSuffix");
# my $RBPfile = "$outputDir/" . "$prefix$siteFileSuffix" . ".all.rbpNum.distribution.txt";
# &getRNase( $RBPfile, $outputDir, $prefix );

# # #15 m6A
# my $extendFile = $outputDir . "/" . "upDownExtend" . $extLen . "nt.bed";
# &RNAmod($extendFile, $outputDir, $prefix);

# # 16 m6A.dastance
my $cuttingSitesBED = "./diffExpSites/" . "$prefix.bed";
# &annotateBedSite( $xmlConfHash, $outputDir, $prefix, $distSiteOutDir );

sub upperCase{
    my ($inputFile) = @_;

    my $outputDir=dirname($inputFile);
    executeCommand("mkdir -p $outputDir") if (!(-e $outputDir));
    open( IN, "<$inputFile" ) || die "can't open the $inputFile\n";

    my $startFile = $outputDir."/". basename($inputFile) . ".upperCase.fa";
    $startFile=~s/\.fa//;
    open( START, ">$startFile" ) || die "can't open the $startFile\n";
    while ( my $line = <IN> ) {
        $line =~ s/\s+$//;
        if ($line ~~/^>/i){
            $line=~s/\(.*\)//;
            print START $line,"\n";
            next;
        }
        print START uc($line),"\n";
    }
    close(IN);
    close(START);
    return ($startFile);
}

sub bedtoolsGetFasta {
    my ( $confHash, $inputFile ) = @_;
    my $outFile = $inputFile;
    $outFile =~ s/\.txt$//;
    $outFile =~ s/\.bed$//;
    $outFile .= ".fa";
    my $commandLine
        = "bedtools2". " getfasta "
        . " -fi " . $confHash->{'genomeFile'}
        . " -bed " . $inputFile
        . " -s -nameOnly "
        . " > ". $outFile;
    &executeCommand("$commandLine");
    return $outFile;
}

sub runShapemapper {
    my ( $samFile, $fa, $outputDir ) = @_;
    my $mutation = $outputDir . "/" . $prefix . ".mutation";
    my $commandLine_1
        = "/public/home/huangjh/softwares/shapemapper-2.1.5/internals/bin/shapemapper_mutation_parser"
        . " -i $samFile -o $mutation";
    &executeCommand($commandLine_1);

    my $count = $outputDir . "/" . $prefix . ".mutation.count";
    my $commandLine_2
        = "/public/home/huangjh/softwares/shapemapper-2.1.5/internals/bin/shapemapper_mutation_counter"
        . " -i $mutation --count_out $count --length 121";
    &executeCommand($commandLine_2);

    my $reactivity = $outputDir . "/" . $prefix . ".mutation.count.reactivity";
    my $commandLine_3
        = "/public/home/huangjh/softwares/shapemapper-2.1.5/internals/bin/make_reactivity_profiles.py"
        . " --fa $fa --counts $count --out $reactivity --rna hsa-napRNA-3";
        # Counted mutations files in the following order: modified, untreated control, denatured control.
    &executeCommand($commandLine_3);

}


sub fetch_bed6 {
    my ( $in ) = @_;
    open (IN, "<$in") or die "Cannot open file: $in $!";
    my $outDir=dirname($in)."/bed6";
    executeCommand("mkdir -p $outDir") if (!(-e $outDir));
    my $out= $in."6";
    open (OUT, ">$out") or die "Cannot open file: $out $!";
    while ( my $line = <IN> ) {
        $line =~ s/\s+$//;
        my @items  = split /\t/, $line;
        my $joinItem1  = join "\t", @items[0..5];
        if ($line =~ /#chrom/){
            print OUT $joinItem1,"\n";
            next;
        } 
        print OUT $joinItem1,"\n";
    }
    close(IN);
    close(OUT);
    return $out;
}

sub bamToSam {
    my ( $bam ) = @_;
    my $bam_sort=$bam.".sorted.bam";
    $bam_sort=~s/\.bam//;
    executeCommand("samtools sort -@ 10 $bam -o $bam_sort &");
    
    my $sam=$bam;
    $sam=~s/\.bam/\.sam/;
    executeCommand("samtools view $bam > $sam");    

    return $sam;    
}

sub checkBamRegion {
    my ( $bam, $bed ) = @_;
    # "bedtools2 intersect -abam /public/home/huangjh/NAPseq/human/cutAdapterBarcodeData/hg38_A13_HepG2_RNAseq_rep1_star_splice/hg38_A13_HepG2_RNAseq_rep1.STAR.Aligned.out.bam -b /public/home/huangjh/NAPseq/human/mergeBamExpression/U15B.bed -ubam > /public/home/huangjh/NAPseq/human/mergeBamExpression/U15B.intersectBam.bed"
    my $intersect_bam=$bed.".intersectBam.bam";
    $intersect_bam=~s/\.bed//g;
    executeCommand("bedtools2 intersect -abam $bam -b $bed -ubam > $intersect_bam");
    my $intersect_sam=$intersect_bam;
    $intersect_sam=~s/\.bam/\.sam/;
    executeCommand("samtools view $intersect_bam > $intersect_sam");
}

sub generateBamIndex {
    my ( $bam ) = @_;
    my $bam_sort=$bam.".sorted.bam";
    $bam_sort=~s/\.bam//;
    executeCommand("samtools sort -@ 10 $bam -o $bam_sort");
    executeCommand("samtools index -@ 10 $bam_sort");
    return $bam_sort;
}

sub plotPofile {
    my ($sortedBamFile, $confHash, $outputDir,  $prefix,
        $readLen,       $platform, $barcodeLen, $primer
    ) = @_;
    my $primerTag = " ";

    my $pairTag = " ";
    if ( $platform eq "pairedEnd" && $sortedBamFile!~/nanopore/i) {
        $pairTag = " -P ";
        $readLen = -1000;
    }
    my $libType = "12";
    if ( $sortedBamFile ~~/_RNAseq_|U15B\.intersectBam/ ) {
        $libType = "21";
    }
    my $primerParameter = "";
    if ( $primer ne "" ) {
        $primerTag
            = " --fa "
            . $confHash->{'genomeFile'}
            . " --fai "
            . $confHash->{'faidxFile'};
        $primerParameter = " --primer " . $primer . " ";
    }
    my $logFile   = $outputDir . "/" . $prefix . ".bamToWig.log.txt";
    my $commandLine
        = $confHash->{'bamToWig'} . " "
        . $confHash->{'bamToWigParameter'}
        . $primerParameter . " " . " -L "
        . $readLen 
        . " --lib-type ". $libType 
        . " "
        . $pairTag . " " . " -u "
        . $barcodeLen . " "
        . $primerTag . " " . " -p "
        . $prefix . " -o "
        . $outputDir
        . " --bam "
        . $sortedBamFile
        . " 2>$logFile";
    &executeCommand($commandLine);
    return 1;
}

sub wigToBigwig {
    my ( $outputDir, $prefix, $confHash ) = @_;
    my $style = '.count';
    if ( $confHash->{'bamToWigParameter'} =~ /rpm/ ) {
        $style = '.rpm';
    }
    my @wiggleFiles = (
        $style . '.minus.end.wg',
        $style . '.minus.coverage.wg',
        $style . '.minus.start.wg',
        $style . '.plus.end.wg',
        $style . '.plus.coverage.wg',
        $style . '.plus.start.wg'
    );
    for ( my $i = 0; $i < scalar(@wiggleFiles); $i++ ) {
        my $wigFile = $outputDir . "/" . $prefix . $wiggleFiles[$i];
        my $bwFile  = $wigFile;
        $bwFile =~ s/\.wg/\.bw/;
        my $commandLine
            = $confHash->{'wigToBigWig'} . " "
            . $wigFile . " "
            . $confHash->{'chromSizeFile'} . " " . " "
            . $bwFile;
        executeCommand($commandLine);
        $commandLine = "rm -f " . $wigFile;
        executeCommand($commandLine);
    }
}


sub runBamExpression {
    my ( $confHash, $annoFile, $bamFile, $outputDir, $prefix ) = @_;
    my $expDir="$outputDir/bamExpression";
    executeCommand("mkdir -p $expDir") if (!(-e $expDir));
    my $outfile   = $expDir . "/" . $prefix . ".bamExpression.txt";
    my $logFile   = $expDir . "/" . $prefix . ".bamExpression.log.txt";
    my $commandLine
        = $confHash->{'bamExpression'} . " "
        . $confHash->{'bamExpressionParameter'} . " "#--pair(paired-end) --noSplice(discard reads span exon-exon) --norm
        . "--anno " . $annoFile
        . " --bam " . $bamFile
        . " -o "    . $outfile
        . " 2>$logFile";
    executeCommand("$commandLine &");
    return $outfile;
}

sub plotMetaGena {
    my ( $confHash, $inSite, $outDir ) = @_;
    my $outMetaGene="$outDir/".basename($inSite).".metagene";
    executeCommand("perl /public/home/huangjh/perlBioTools/myPerlScript/metagene/bedBinDistribution.pl --input $inSite -bed6 /public/home/huangjh/perlBioTools/myPerlScript/metagene/gencode.v30.chr_patch_hapl_scaff.bed6 -o $outMetaGene");
}

sub annotateSitesToBed {
    my ( $confHash, $annotateSiteFile, $outDir, $minVal) = @_;
    my $siteDir = $confHash->{'splicingDir_V30'};
    executeCommand("mkdir -p $outDir") if (!(-e $outDir));
    my $prefix=basename($annotateSiteFile);
    my $knownSite = "$siteDir/" . (split/\./,$prefix)[1] . "_sites.bed";
    my ( $exonPtr, $genePtr, $linePtr ) = &fetchKnownPos($knownSite);
    my ( $detectSite, $unDetectSite )   = &fetchAnnotateSites( $annotateSiteFile, $exonPtr, $genePtr, $linePtr, $prefix, $outDir, $minVal);
    return ( $detectSite, $unDetectSite ); 
}

sub fetchKnownPos {
    my $knownSite = shift @_;
    my %lineHash = ();
    my %exonHash = ();
    my %geneHash = ();
    open( EXON, "<$knownSite" ) || die "can't open the $knownSite\n";
    while ( my $line = <EXON> ) {
        $line =~ s/\s+$//;
        my @items  = split /\t/, $line;
        my $chrom  = $items[0];
        my $start  = $items[1];
        my $end    = $items[2];
        my $name   = $items[3];
        my $score  = $items[4];
        my $strand = $items[5];
        my $position    = $chrom . ":" . $start . ":" . $end . ":" . $strand;
        $exonHash{$name} = $position;
        $lineHash{$name} = $line;
        $geneHash{$name}{$position}++;
    }
    close(EXON);
    return ( \%exonHash, \%geneHash, \%lineHash );
}

sub fetchAnnotateSites {
    my ( $annotateSitesFile, $exonPtr, $genePtr,$linePtr, $prefix, $outDir, $minVal ) = @_;
    my $colIdx      = 3;
    my %annotateSitesNameHash = ();
    my @headItems   = ();
    open( NAP, "<$annotateSitesFile" ) || die "can't open the $annotateSitesFile\n";
    while ( my $line = <NAP> ) {
        $line =~ s/\s+$//;
        my @items = split /\t/, $line;
        if ( $line =~ /^DistanceFromCenter/ ) {
            for ( my $i = 0; $i < scalar(@items); $i++ ) {
                push @headItems, $items[$i];
            }
            next;
        }
        if ( $line =~ /^0\t\d*/ ) {
            for ( my $i = $colIdx; $i < scalar(@items); $i++ ) {
                $annotateSitesNameHash{ $headItems[$i] } = $items[$i];####score
            }
            next;
        }
    }
    close(NAP);

    my $i = 1;
    my $sbEndFile= "$outDir/" . $prefix . ".detectSite.bed";
    open( SBENDSITE, ">$sbEndFile" );
    print SBENDSITE "#chrom\tchromStart\tchromEnd\tname\tendRPM\tstrand","\n";
    my %geneHash       = ();
    my %stableEndPosHash     = ();
    my %allPosHash     = ();
    foreach my $stableName ( sort { $annotateSitesNameHash{$b} <=> $annotateSitesNameHash{$a} }keys %annotateSitesNameHash ) {
        $allPosHash{ $exonPtr->{$stableName} }++;
        if ( $annotateSitesNameHash{$stableName} >= $minVal && defined( $exonPtr->{$stableName} ) ) {# our $minVal = 0.5;
            my $geneName = $stableName;
            my $position = $exonPtr->{$stableName};
            $stableEndPosHash{$geneName} = $position;
            $geneHash{ $stableName }++;
            my @endInfo = split /\:/, $exonPtr->{$stableName};
            my $chrom   = $endInfo[0];
            my $name    = $stableName;
            my $strand  = $endInfo[3];
            my $knownSiteLine = $linePtr->{$name};
            print SBENDSITE $chrom, "\t", $endInfo[1], "\t", $endInfo[2], "\t", $stableName,
                "\t$annotateSitesNameHash{$stableName}\t", $strand,"\n";
            $i++;
        }
    }
    close(SBENDSITE);

    my $other= "$outDir/" . $prefix . ".unDetectSite.bed";
    open( OTHER, ">$other" );
    $i = 1;
    foreach my $name ( keys %{$exonPtr} ) {
        if ( !defined( $stableEndPosHash{$name} ) ) {
            my $knownSiteLine = $linePtr->{$name};
            print OTHER $knownSiteLine, "\n";
        }
    }
    close(OTHER);
    return ( $sbEndFile, $other );
}

sub annotateChimeras {
    my ( $confHash, $targetFile, $outputDir, $prefix ) = @_;
    $outputDir .= "/";
    my $targetBed  = $targetFile . ".candidate.bed";
    my %targetHash = ();
    open( TARGET, "<$targetFile" ) || die "can't open the $targetFile\n";
    open( OUTPUT, ">$targetBed" ) || die "can't open the $targetBed\n";
    while ( my $line = <TARGET> ) {
        $line =~ s/\s+$//;
        next if ( $line =~ /lChrom/i && $line =~ /rChrom/i );
        my @items = split /\t/, $line;
        my $keyVal = $items[3] . "split" . $items[9];
        $targetHash{$keyVal} = $line;
        print OUTPUT $items[0], "\t", $items[1], "\t", $items[2], "\t", $items[3], "\t", $items[4], "\t", $items[5], "\n";
        print OUTPUT $items[6], "\t", $items[7], "\t", $items[8], "\t", $items[9], "\t", $items[10], "\t", $items[11], "\n";
    }
    close(TARGET);
    close(OUTPUT);

    ######嵌合体分开两端分别注释
    my $conMap = &getResultConVals_rriScan( $targetBed, $confHash );
    my $annotatedMap  = &getAnnotateResults( $targetBed, $confHash );
    my %annotatedHash = %{$annotatedMap};
    my %conHash       = %{$conMap};

    my $annoChimeraFile = $outputDir . $prefix . ".rriScan.bedAnnotator.chimera.txt";
    open( OUTPUT, ">$annoChimeraFile" ) || die "can't open $annoChimeraFile\n";
    foreach my $name ( sort keys %targetHash ) {
        my @names = split /split/, $name;
        print OUTPUT $targetHash{$name}, "\t", $conHash{ $names[0] }, "\t", $conHash{ $names[1] }, "\t", $annotatedHash{ $names[0] }, "\t", $annotatedHash{ $names[1] }, "\n";
    }
    close(OUTPUT);
    my $commandLine = "rm -f " . $targetBed;
    executeCommand($commandLine);
    return $annoChimeraFile;
}


sub filterCandidates {
    my ( $confHash, $annoFile, $outputDir, $prefix ) = @_;
    my %rriTypeHash = ();
    my %rriAnnoHash = ();
    my $allNum      = 0;
    my $filterNum   = 0;

    my $filterFile
        = $outputDir . "/" . $prefix . ".rriScan.chimara.filter.txt";
    open( OUTPUT, ">$filterFile" ) || die "can't open the $filterFile\n";
    open( ANNO,   "<$annoFile" )   || die "can't open the $annoFile\n";

    my $header
        = "lChrom\tlChromStart\tlChromEnd\tlName\tlScore\tlStrand\trChrom\trChromStart\trChromEnd\trName\trScore\trStrand\t";
    $header
        .= "lociNum\tgapDist\treadSeq\tchimeraSeq\tchimeraStruct\tMFE\trriType\t";
    $header
        .= "lAlignSeq\tpairs\trAlignSeq\tpairNum\talignScore\tloReadNum\troReadNum\tlConservationScore\trConservationScore\t";
    $header
        .= "lAnnoInfo\tlAnnoGene\tlAnnoType\trAnnoInfo\trAnnoGene\trAnnoType";

    print OUTPUT $header, "\n";

    while ( my $line = <ANNO> ) {
        $line =~ s/\s+$//;
        my @items    = split /\t/, $line;
        my $mfe      = $items[17];
        my $pairNum  = $items[22];
        my $alnScore = $items[23];
        my $lChrom   = $items[0];
        my $rChrom   = $items[6];
        my $rriType  = $items[18];
        my $rType    = $items[-1];
        my $lType    = $items[-4];
        $allNum++;

        if (   $pairNum >= $minPair
            && $mfe < $maxMfe
            && $alnScore >= $minScore )
        {
            $filterNum++;
            $rriTypeHash{$rriType}++;
            my $annoKey
                = $lType le $rType
                ? $lType . "#" . $rType
                : $rType . "#" . $lType;
            $rriAnnoHash{$annoKey}++;
            print OUTPUT $line, "\n";
        }
    }
    close(ANNO);
    close(OUTPUT);

    my $rriTypeFile
        = $outputDir . "/" . $prefix . ".rriScan.chimera.rriType.txt";
    open( RRITYPE, ">$rriTypeFile" ) || die "can't open the $rriTypeFile\n";
    print RRITYPE "type\tnumber\n";
    foreach my $key (
        sort { $rriTypeHash{$b} <=> $rriTypeHash{$a} }
        keys %rriTypeHash
        )
    {
        print RRITYPE $key, "\t", $rriTypeHash{$key}, "\n";
    }
    close(RRITYPE);

    my $annoTypeFile
        = $outputDir . "/" . $prefix . ".rriScan.chimera.annoType.txt";
    open( ANNOTYPE, ">$annoTypeFile" )
        || die "can't open the $annoTypeFile\n";
    print ANNOTYPE "type\tnumber\n";
    foreach my $key (
        sort { $rriAnnoHash{$b} <=> $rriAnnoHash{$a} }
        keys %rriAnnoHash
        )
    {
        print ANNOTYPE $key, "\t", $rriAnnoHash{$key}, "\n";
    }
    close(ANNOTYPE);

    return $filterNum;
}

sub getChimericRatio {
    my ( $confHash, $filterNum, $outputDir, $prefix ) = @_;

    my $mapReadNum  = 0;
    my $uniNum      = 0;
    my $mulNum      = 0;
    my $chimericNum = 0;
    my $gapNum      = 0;

    my $starSummaryFile = $outputDir . "/" . $prefix . ".STAR.Log.final.out";
    open( MAP, "<$starSummaryFile" )
        || die "can't open the $starSummaryFile\n";
    while ( my $line = <MAP> ) {
        $line =~ s/\s+$//;
        if ( $line =~ /Uniquely\s+mapped\s+reads\s+number\s+\|\s+(\d+)/ ) {
            $uniNum = $1;
        }
        if ( $line
            =~ /Number\s+of\s+reads\s+mapped\s+to\s+multiple\s+loci\s+\|\s+(\d+)/
            )
        {
            $mulNum = $1;
        }
        if ( $line =~ /Number\s+of\s+chimeric\s+reads\s+\|\s+(\d+)/ ) {
            $chimericNum = $1;
        }
    }
    close(MAP);

    my $rriScanLogFile = $outputDir . "/" . $prefix . ".junction.rriScan.log";
    open( RRISCAN, "<$rriScanLogFile" )
        || die "can't open the $rriScanLogFile\n";
    while ( my $line = <RRISCAN> ) {
        $line =~ s/\s+$//;
        if ( $line =~ /identify\s+(\d+)\s+chimeric.*bam/ ) {
            $gapNum = $1;
        }
    }
    close(RRISCAN);

    my $ratioFile = $outputDir . "/" . $prefix . ".rriScan.chimericRatio.txt";
    open( RATIO, ">$ratioFile" ) || die "can't open the $ratioFile\n";
    print RATIO
        "filterNum\tchimericNum\tgapNum\tuniqueNum\tmultipleNum\tfilterPercent\tallPercent\n";
    print RATIO $filterNum, "\t", $chimericNum, "\t", $gapNum, "\t", $uniNum,
        "\t", $mulNum, "\t",
        $filterNum / ( $filterNum + $uniNum + $mulNum ) * 100,
        "\t",
        ( $gapNum + $chimericNum )
        / ( $gapNum + $chimericNum + $uniNum + $mulNum )
        * 100, "\n";
    close(RATIO);

    return $ratioFile;
}


sub mergeIntron {
    my ( $confHash, $ncrnaFile, $outputDir, $prefix ) = @_;
    my %mergeHash = ();
    my $header    = "";
    open( NCRNA, "<$ncrnaFile" ) || die "can't open the $ncrnaFile\n";

    #chr1:15560083-15560125(+)
    while ( my $line = <NCRNA> ) {
        $line =~ s/\s+$//;
        if ( $line =~ /chrom/ ) {
            $header = $line;
            next;
        }
        my @items      = split /\t/, $line;
        my $chrom      = $items[0];
        my $chromStart = $items[1];
        my $chromEnd   = $items[2];
        my $name       = $items[3];
        my $readNum    = $items[4];
        my $strand     = $items[5];
        my $starts     = $items[11];
        my $keyStr
            = $chrom . ":"
            . $chromStart . ":"
            . $chromEnd . ":"
            . $strand . ":"
            . $starts;

        if ( defined( $mergeHash{$keyStr} ) ) {
            $mergeHash{$keyStr}->[0] += 1;
        }
        else {
            $mergeHash{$keyStr} = [ 1, $line ];
        }
    }    # while end
    close(NCRNA);

    my $mergeFile = $outputDir . "/" . $prefix . ".intronScan.mergeIntron.txt";
    open( MERGE, ">$mergeFile" ) || die "can't open the $mergeFile\n";
    print MERGE $header, "\n";
    foreach my $key ( sort keys %mergeHash ) {
        my @items = split /\t/, $mergeHash{$key}->[1];
        $items[4] = $mergeHash{$key}->[0];
        print MERGE join( "\t", @items ), "\n";
    }
    close(MERGE);
    return $mergeFile;
}

sub alignReadToGenome_rriScan {
    my ($read1File, $read2File, $confHash, $outputDir,
        $prefix,    $aligner,   $platform, $seqFormat
    ) = @_;

    my $readFiles        = $read1File;
    my $starOutputPrefix = $outputDir . "/" . $prefix . ".STAR.";
    if ( $confHash->{'starParameter_rriScan'} =~ /Local/ ) {
        $starOutputPrefix = $outputDir . "/" . $prefix . ".Local.STAR.";
    }
    if ( $platform eq "pairedEnd" ) {
        $readFiles = "" . $read1File . " " . $read2File . " ";
    }

    my $commandLine = &alignReadStar(
        $confHash->{'starProgram'},
        $confHash->{'starParameter_rriScan'},
        $confHash->{'starIdx'},
        $readFiles, $starOutputPrefix
    );
    &executeCommand($commandLine);

    my $sortBamFile = $starOutputPrefix . "Aligned.sortedByCoord.out.bam";
    my $mate1File   = $starOutputPrefix . "Unmapped.out.mate1";
    my $mate2File   = $starOutputPrefix . "Unmapped.out.mate2";

    if ( -e $mate1File ) {
        my $gzipCommandLine = "nohup gzip -f $mate1File &";
        &executeCommand($gzipCommandLine);
    }
    if ( -e $mate2File ) {
        my $gzipCommandLine = "nohup gzip -f $mate2File &";
        &executeCommand($gzipCommandLine);
    }
    $mate1File .= ".gz";
    $mate2File .= ".gz";

    return ( $sortBamFile, $mate1File, $mate2File );
}

sub runRriScan {
    my ( $confHash, $sortedBamFile, $outputDir, $prefix, $junFile ) = @_;
    my $ncapSeeker = $outputDir . "/" . $prefix . ".rriScan";
    my $commandLine
        = $confHash->{'rriScan'} . " "#/public/home/jianhua/bioTools/bin/ncapSeeker
        . $confHash->{'rriScanParameter'} . " "#-t 5 -c 2 -F 2 -e 5 -s --norm --pair
        . " --fa "
        . $confHash->{'genomeFile'}#/public/home/jianhua/genome/human/hg38/WholeGenomeFasta/hg38.fa
        . " --fai "
        . $confHash->{'faidxFile'}#/public/home/jianhua/genome/human/hg38/WholeGenomeFasta/hg38.fa.fai
        . " --bam "
        . $sortedBamFile . " "
        . " --jun "
        . $junFile . " "
        . " >$ncapSeeker";
    &executeCommand($commandLine);
    return $ncapSeeker;
}

sub runIntronScan {
    my ( $confHash, $sortedBamFile, $outputDir, $prefix ) = @_;
    my $ncapSeeker = $outputDir . "/" . $prefix . ".intronScan";
    my $commandLine
        = $confHash->{'intronScan'} . " "#/public/home/jianhua/bioTools/bin/ncapSeeker
        . $confHash->{'intronScanParameter'} . " "#-t 5 -c 2 -F 2 -e 5 -s --norm --pair
        . " --fa "
        . $confHash->{'genomeFile'}#/public/home/jianhua/genome/human/hg38/WholeGenomeFasta/hg38.fa
        . " --fai "
        . $confHash->{'faidxFile'}#/public/home/jianhua/genome/human/hg38/WholeGenomeFasta/hg38.fa.fai
        . " --bam "
        . $sortedBamFile . " "
        . " >$ncapSeeker";
    &executeCommand($commandLine);
    return $ncapSeeker;
}


sub RNAmod{
    my ( $extendFile, $outputDir, $prefix ) = @_;
    my $m6AFile = "/public/home/huangjh/genome/modSites/human_all_RNA_modification_sites.bed6";
    my $intersectFile = "$outputDir/" . $prefix . "intersectRNAmod.bed";
    my $statisticFile = "$outputDir/" . $prefix . "intersectRNAmod.statistic.bed";
    executeCommand("intersectBed -a $extendFile -b $m6AFile -wa -wb > $intersectFile");

    my %modHash = ();
    my %napID = ();
     
    open( ANNO, "<$intersectFile" );

    while ( my $line = <ANNO> ) {

        $line =~ s/\s+$//;
        next if ( $line =~ /#/ );

        my @lineItem = split /\t/, $line;
        
        my $napID = $lineItem[3];

        my $RNAmod = $lineItem[9];

        $modHash{ $napID } .= $RNAmod . ",";
    }

    open(OUT, ">$statisticFile");

    print OUT "napID\t", "RNAmod\t", "m1A\t", "m5C\t", "m6A\t", "m7G\t", "Nm\t", "Pseudo\t", "RNAediting\t", "otherMod\n";

    foreach my $napID ( sort keys %modHash ){

        my $m1ACounts = ($modHash{ $napID } =~ s/m1A_site/m1A_site/g) + 0;
        my $m5CCounts = ($modHash{ $napID } =~ s/m5C_site/m5C_site/g) + 0;
        my $m6ACounts = ($modHash{ $napID } =~ s/m6A_site/m6A_site/g) + 0;
        my $m7GCounts = ($modHash{ $napID } =~ s/m7G_site/m7G_site/g) + 0;
        my $NmCounts  = ($modHash{ $napID } =~ s/Nm_site/Nm_site/g) + 0;
        my $PseudoCounts = ($modHash{ $napID } =~ s/Pseudo_site/Pseudo_site/g) + 0;
        my $editingCounts = ($modHash{ $napID } =~ s/RNA-editing_site/RNA-editing_site/g) + 0;
        my $allCounts = ($modHash{ $napID } =~ s/,/,/g) + 0;

        my $otherModCounts = $allCounts-$m1ACounts-$m5CCounts-$m6ACounts-$m7GCounts-$NmCounts-$PseudoCounts-$editingCounts;

        print OUT "$napID\t", "$modHash{ $napID }\t", "$m1ACounts\t", "$m5CCounts\t", "$m6ACounts\t", "$m7GCounts\t", "$NmCounts\t", "$PseudoCounts\t", "$editingCounts\t", "$otherModCounts\n";
    }

    close(OUT);
    close(ANNO);
}

sub getRNase{
    my ( $RBPfile, $outputDir, $prefix ) = @_;#all.rbpNum.distribution.txt

    my %RBPhash = ();
    my %endoHash = ();
    my %exoHash = ();

    open( RBP, "<$RBPfile" ) || die "can't open the $RBPfile\n";
    while ( my $line = <RBP> ) {

        $line =~ s/\s+$//;
        my @items = split /\t/, $line;
        my $RBP = $items[1];
        my $RBPbindingPercent = $items[3];

        if ($line =~ /intronPercent/i) {
             next;
         }  
        
        else{
            $RBPhash{$RBP} = $RBPbindingPercent;
        }

    }

    open( UNIPROT, "</public/home/huangjh/genome/UniProtKB/human_ALLproteinFunction_UniProtKB.txt" );

    my $outFile = "$outputDir/" . $prefix . ".RNase.txt";
    open( OUTPUT1, ">$outFile" );



    while ( my $line = <UNIPROT> ) {
        $line =~ s/\s+$//;
        my @items = split /\t/, $line;
        my $geneName = $items[5];
        my $function = $items[7];

        if ($line =~ /#stringID/i) {
             print OUTPUT1 "#RNase\tRBP\tRBPbindingPercent", "\t", $line, "\n";
        }

        else{
            foreach my $rbp (keys %RBPhash){

                #endonuclease;endoribonuclease
                #exonuclease;exoribonuclease
                #RNase;Ribonuclease  #cleavage
                if ($geneName =~ m/$rbp/i) {
                    if ($function =~ /endonuclease/i || $function =~ /endoribonuclease/i){
                        $endoHash{$geneName} = "endonuclease" . "\t" . $rbp . "\t" . $RBPhash{$rbp} . "\t" . $line;
                    }

                    if ($function =~ /exonuclease/i || $function =~ /exoribonuclease/i){
                        $exoHash{$geneName} = "exonuclease" . "\t" . $rbp . "\t" . $RBPhash{$rbp} . "\t" . $line;
                    }


                }
            }                

         }
    }

    foreach my $RNase (keys %endoHash){
        print OUTPUT1 $endoHash{$RNase}, "\n"; 
    }
    
    foreach my $RNase (keys %exoHash){
        print OUTPUT1 $exoHash{$RNase}, "\n"; 
    }



    close(UNIPROT);
    close(OUTPUT1);
}




#fetch cuttingSite

sub cuttingSiteFilter{
    my ( $mergeBed ) = @_;
    my $cuttingSiteFile = "$outputDir/$prefix.STAR.Aligned.out.bam.startEndSeeker.CutSitesAnnotator.rmAnnoCutSites.merge.filterCutSites.bed";

    open( INPUT, "<$mergeBed" );
    open( OUTPUT, ">$cuttingSiteFile" );

        while ( my $line = <INPUT> ) {
            $line =~ s/\s+$//;
            my @items = split /\s+/, $line;
            my $mergesites      = $items[3];
            my $mergeNums       = $items[4];

            if ($mergesites =~ /start/i && $mergesites =~ /end/i){
                my @sites  = split /;/, $mergesites;
                if (  $sites[-1] =~ /start/){
                    print OUTPUT "$line\n";   
                }
            }
        }

    close(INPUT);
    close(OUTPUT);

    return $cuttingSiteFile
}


#1
sub computeStartEndTag {
    my ( $confHash, $readFile, $outputDir, $prefix ) = @_;
    my $outFile = $outputDir . "/" . $prefix . ".startEndTagSummary.txt";#/public/home/huangjh/NAPseq/human/NCAP_seq_condition_data/cutAdapterBarcodeData/hg38_AN-1201213_ADR_NAP-seq-rep1_star_splice
    open( OUTPUT, ">$outFile" ) || die "can't open the $outFile\n";
    print OUTPUT "#sampleName", "\t", "totalRead", "\t", "p3Read", "\t",
        "p5Read",
        "\t", "p8Read", "\t", "othRead", "\t", "startEnds", "\t",
        "percentage", "\n";
    warn "# deal with the file: " . $readFile, "\n";
    my $totalRead = 0;
    my $p3Read    = 0;
    my $p5Read    = 0;
    my $p8Read    = 0;
    my $othRead   = 0;
    open( my $rfh, "gzip -dc \"$readFile\" |" )#-c将输出写到标准输出并保留原有文件;-d将压缩文件解压。
        || die("Could not open gzipped read file: $readFile \n");

    while ( my $seqName = <$rfh> ) {
        my $seq     = <$rfh>;
        my $qltName = <$rfh>;
        my $qlt     = <$rfh>;
        $seqName =~ s/\s+$//;
        $totalRead++;
        if ( $seqName =~ /\:3$/ ) {
            $p3Read++;
        }
        elsif ( $seqName =~ /\:5$/ ) {
            $p5Read++;
        }
        elsif ( $seqName =~ /\:8$/ ) {
            $p8Read++;
        }
        else {
            $othRead++;
        }
    }
    close($rfh);
    print OUTPUT $readFile, "\t", $totalRead, "\t", $p3Read, "\t",
        $p5Read, "\t",
        $p8Read, "\t", $othRead, "\t", ( $p3Read + $p5Read + $p8Read ),
        "\t", ( $p3Read + $p5Read + $p8Read ) / $totalRead * 100,
        "\n";

    close(OUTPUT);

    return $outFile;
}


#2
sub alignReadToNapRNAs {
    my ($read1File, $read2File, $confHash, $outputDir,
        $prefix,    $aligner,   $platform, $seqFormat
    ) = @_;

    my $readFiles        = $read1File;
    my $starOutputPrefix = $outputDir . "/" . $prefix . ".STAR.";
    # if ( $confHash->{'starParameter'} =~ /Local/ ) {
    #     $starOutputPrefix = $outputDir . "/" . $prefix . ".Local.STAR.";
    # }
    if ( $platform eq "pairedEnd" ) {
        $readFiles = "" . $read1File . " " . $read2File . " ";
    }

    my $commandLine = &alignReadStar(
        $confHash->{'starProgram'},
        $confHash->{'starParameter_shapemapper'},
        $confHash->{'starIdx_napRNAs'},
        $readFiles, $starOutputPrefix
    );
    &executeCommand($commandLine);

    my $sortBamFile = $starOutputPrefix . "Aligned.out.bam";
    my $mate1File   = $starOutputPrefix . "Unmapped.out.mate1";
    my $mate2File   = $starOutputPrefix . "Unmapped.out.mate2";

    if ( -e $mate1File ) {
        my $gzipCommandLine = "nohup gzip $mate1File &";
        &executeCommand($gzipCommandLine);
    }
    if ( -e $mate2File ) {
        my $gzipCommandLine = "nohup gzip $mate2File &";
        &executeCommand($gzipCommandLine);
    }
    $mate1File .= ".gz";
    $mate2File .= ".gz";

    return ( $sortBamFile, $mate1File, $mate2File );
}

#2
sub alignReadToGenome {
    my ($read1File, $read2File, $confHash, $outputDir,
        $prefix,    $aligner,   $platform, $seqFormat
    ) = @_;

    my $readFiles        = $read1File;
    my $starOutputPrefix = $outputDir . "/" . $prefix . ".STAR.";
    # if ( $outputDir =~ /SHAPE/i ) {
    #     $starOutputPrefix = $outputDir . "/" . $prefix . "STAR.genome.";
    # }
    if ( $platform eq "pairedEnd" ) {
        $readFiles = "" . $read1File . " " . $read2File . " ";
    }

    my $commandLine = &alignReadStar(
        $confHash->{'starProgram'},
        $confHash->{'starParameter'},
        $confHash->{'starIdx'},
        $readFiles, $starOutputPrefix
    );
    &executeCommand($commandLine);

    my $sortBamFile = $starOutputPrefix . "Aligned.out.bam";
    my $mate1File   = $starOutputPrefix . "Unmapped.out.mate1";
    my $mate2File   = $starOutputPrefix . "Unmapped.out.mate2";

    if ( -e $mate1File ) {
        my $gzipCommandLine = "nohup gzip $mate1File &";
        &executeCommand($gzipCommandLine);
    }
    if ( -e $mate2File ) {
        my $gzipCommandLine = "nohup gzip $mate2File &";
        &executeCommand($gzipCommandLine);
    }
    $mate1File .= ".gz";
    $mate2File .= ".gz";

    return ( $sortBamFile, $mate1File, $mate2File );
}

#3
sub runNcapSeeker {
    my ( $confHash, $sortedBamFile, $outputDir, $prefix ) = @_;
    my $ncapSeeker = $outputDir . "/" . $prefix . ".ncapSeeker";
    my $commandLine
        = $confHash->{'ncapSeeker'} . " "#/public/home/jianhua/bioTools/bin/ncapSeeker
        . $confHash->{'ncapSeekerParameter'} . " "#-t 5 -c 2 -F 2 -e 5 -s --norm --pair
        . " --fa "
        . $confHash->{'genomeFile'}#/public/home/jianhua/genome/human/hg38/WholeGenomeFasta/hg38.fa
        . " --fai "
        . $confHash->{'faidxFile'}#/public/home/jianhua/genome/human/hg38/WholeGenomeFasta/hg38.fa.fai
        . " --input "
        . $sortedBamFile . " "
        . " >$ncapSeeker";
    if ($sortedBamFile~~/nanopore/i){
        $commandLine
        = $confHash->{'ncapSeeker'} . " "
        . $confHash->{'ncapSeekerParameter_nanopore'} . " "
        . " --fa "
        . $confHash->{'genomeFile'}
        . " --fai "
        . $confHash->{'faidxFile'}
        . " --input "
        . $sortedBamFile . " "
        . " >$ncapSeeker";
    }
    &executeCommand($commandLine);
    return $ncapSeeker;
}

sub getResultConVals_rriScan {
    my ( $bedFile, $confHash ) = @_;
    my $conValFile = &getEvoConVals(
        $confHash->{'bedConservation'},
        $confHash->{'bedConservationParameter'},
        $confHash->{'conservationProfile'}, $bedFile
    );
    my %conHash = ();
    open( ANNO, "<$conValFile" )
        || die "can't open $conValFile\n";
    while ( my $line = <ANNO> ) {
        $line =~ s/\s+$//;
        my @items = split /\t/, $line;
        $conHash{ $items[3] } = $items[$#items];
    }
    close(ANNO);
    my $commandLine = "rm -f $conValFile";
    &executeCommand($commandLine);
    return \%conHash;
}

sub getAnnotateResults {
    my ( $bedFile, $confHash ) = @_;
    my $annotatedFile = $bedFile . ".bedAnnotator.txt";
    my $commandLine
        = $confHash->{'bedAnnotator'} . " "
        . $confHash->{'bedAnnotatorParameter'}
        . " --anno "
        . $confHash->{'exonAnnotation'}
        . " --bed "
        . $bedFile . " > "
        . $annotatedFile;
    &executeCommand($commandLine);

    my %annotatedHash = ();
    open( ANNO, "<$annotatedFile" )
        || die "can't open $annotatedFile\n";
    while ( my $line = <ANNO> ) {
        $line =~ s/\s+$//;
        my @items      = split /\t/, $line;
        my $geneName   = $items[-7];
        my $geneType   = $items[-6];
        my @geneInfo   = split /\|/, $geneName;
        my $geneSymbol = $geneName;
        if ( scalar(@geneInfo) > 3 ) {
            $geneSymbol = $geneInfo[3];
        }
        $annotatedHash{ $items[3] }
            = $geneName . "\t" . $geneSymbol . "\t" . $geneType;
    }
    close(ANNO);
    $commandLine = "rm -f $annotatedFile";
    &executeCommand($commandLine);
    return \%annotatedHash;
}

#4
sub getResultConVals {
    my ( $confHash, $inputFile ) = @_;
    my $conValFile = &getEvoConVals(#？？？这子程序没见到有定义
        $confHash->{'bedConservation'},#/public/home/jianhua/bioTools/bin/bedConservation，bedConservation什么程序？
        $confHash->{'bedConservationParameter'},#>-c -1000
        $confHash->{'conservationProfile'}, $inputFile#/public/home/jianhua/genome/human/hg38/conservation/hg38.phyloP100way.bw
    );
    return $conValFile;
}

#5.1 普通的bedannotator
sub annotateResults {
    my ( $confHash, $inputFile ) = @_;
    my $annotatedFile = $inputFile . ".bedAnnotator";

    my $commandLine
        = $confHash->{'bedAnnotator'} . " "##/public/home/jianhua/bioTools/bin/bedAnnotator
        . $confHash->{'bedAnnotatorParameter'}##-s 1 #1 is same strand
        . " --anno "
        . $confHash->{'exonAnnotation'}##/public/home/jianhua/genome/human/hg38/annotations/gencodeV30/annotationBed/hg38.genecode.v30.tRNA.snoRNA.miRNA.rmsk.chrM.RefseqNcrna.curated.exonFeatures.bed6
        . " --bed "
        . $inputFile . " > "
        . $annotatedFile;
    &executeCommand($commandLine);
    return $annotatedFile;
}


#5.2 普通的bedannotator总结
sub annotateResultsSummary {
    my ( $confHash, $inputFile ) = @_;
    my $annotatedFileSummary = $inputFile . ".bedAnnotator.summary.txt";


    my $commandLine
        = "/public/home/huangjh/C/bedAnnotatorTools/bin/bedAnnotator_RNAseH" . " "
        . $confHash->{'bedAnnotatorParameter'}##-s 1 #1 is same strand
        . " --summary"
        . " --anno "
        . $confHash->{'exonAnnotation_revisedSnoRA37AB'}#注释文件，按需改
        . " --bam "
        . $inputFile . " > "
        . $annotatedFileSummary;

    &executeCommand($commandLine);
    return $annotatedFileSummary;
}

sub pieplot_common {
    my ( $confHash, $summaryFile, $outDir, $annoIdx, $plotTopNums ) = @_;
    my $sum=0;
    executeCommand("mkdir $outDir/pieplot") if (!(-e "$outDir/pieplot"));
    my $outFile= "$outDir/pieplot/".basename($summaryFile). ".pieplot". "_top"."$plotTopNums".".txt";

    my %typeHash = ();
    open( TYPE, "<$summaryFile") || die "can't open $summaryFile";
    open( OUT, ">$outFile") || die "can't open $outFile";
    print OUT "#lable\t", "percent\n";

    while ( my $line = <TYPE> ){
        $line=~s/\s+$//;
        next if ( $line =~ /\#/ );
        my @items = split /\t/, $line;

        my $annoType = $items[$annoIdx];
        my $readNum=$items[1];

        if( $annoType=~/intron/i ){
            $typeHash{"intron"}+=$readNum;
            $sum+=$readNum; 
        }

        elsif( $annoType=~/CDS/i ){
            $typeHash{"CDS"}+=$readNum;
            $sum+=$readNum; 
        }
        
        elsif( $annoType=~/repeat/i ){
            $typeHash{"repetitive-element"}+=$readNum;
            $sum+=$readNum;
        }

        elsif( $annoType=~/rRNA/i ){
            $typeHash{"rRNA"}+=$readNum;
            $sum+=$readNum; 
        }

        elsif( $annoType=~/intergenic/i ){
            $typeHash{"intergenic"}+=$readNum;
            $sum+=$readNum; 
        }

        elsif( $annoType=~/RMRP/i ){
            $typeHash{"RMRP"}+=$readNum;
            $sum+=$readNum; 
        }
        elsif( $annoType=~/RPPH1/i ){
            $typeHash{"RPPH1"}+=$readNum;
            $sum+=$readNum; 
        }
        elsif( $annoType=~/RPPH1/i ){
            $typeHash{"RPPH1"}+=$readNum;
            $sum+=$readNum; 
        }
        elsif( $annoType=~/SNORA73A/i ){
            $typeHash{"SNORA73A"}+=$readNum;
            $sum+=$readNum;
        }
        elsif( $annoType=~/SNORA73B/i ){
            $typeHash{"SNORA73B"}+=$readNum;
            $sum+=$readNum;
        }
        elsif( $annoType =~ /snoRNA/i
                || $annoType =~ /CDBox/i
                || $annoType =~ /SNORD/i
                || $annoType =~ /scaRNA/i
                || $annoType =~ /SNORA/i
                || $annoType =~ /HAcaBox/i
                || $annoType =~ /snRNA/i
                || $annoType =~ /miRNA/i
                || $annoType =~ /scRNA/i
                || $annoType =~ /misc_RNA/i
                || $annoType =~ /mascRNA-menRNA/i
                || $annoType =~ /ribozyme/i
                || $annoType =~ /sRNA/i 
                || $annoType =~ /sncRNA/i
                || $annoType =~ /lncRNA/i
                || $annoType =~ /Pol3_7SL_srpRNA/i){
                $typeHash{"ncRNA"}+=$readNum;
                $sum+=$readNum;
        }

        elsif( $annoType=~/UTR3/i ){
            $typeHash{"UTR3"}+=$readNum;
            $sum+=$readNum;
        }

        elsif( $annoType=~/UTR5/i ){
            $typeHash{"UTR5"}+=$readNum;
            $sum+=$readNum;
        }

        else{
            $typeHash{"others"}+=$readNum;
            # $typeHash{"$annoType"}+=$readNum;
            $sum+=$readNum;
        }
    }

    # $typeHash{"1tRNA"}=0 if (!(defined $typeHash{"1tRNA"}));
    my $calculator =0;
    my @plot_arry=();
    foreach my $key (sort { $typeHash{$b}<=> $typeHash{$a} } keys %typeHash){
        if ($calculator<$plotTopNums){
            my $percent = sprintf "%.30f", $typeHash{$key}/$sum;
            print OUT "$key\t", "$typeHash{$key}\n";
            push @plot_arry,$key;
        }
        $calculator++;
    }
    
    close(TYPE);
    close(OUT);
    undef %typeHash;
    executeCommand("Rscript /public/home/huangjh/perlBioTools/ncapSeeker/pie_plot.R -f $outFile -o $outFile");

}

#pieplot
sub mySummaryToPlotPie {
    my ( $annotatedFileSummary, $prefix ) = @_;
    my %typeHash = ();
    open( TYPE, "<$annotatedFileSummary") || die "can't open $annotatedFileSummary";

    while ( my $line = <TYPE> ){
            next if ( $line =~ /\#/ );
            my @items = split /\t/, $line;

            my $annoType = $items[0];
            my $percentage = $items[2];

            if ($annoType =~ /utr3/){
                $typeHash{"utr3"} += $percentage;
            }

            elsif ($annoType =~ /utr5/){
                $typeHash{"utr5"} += $percentage;
            }

            elsif ($annoType =~ /repeatMask/){
                $typeHash{"repeatMask"} += $percentage;
            }

            elsif ($annoType =~ /cds/){
                $typeHash{"cds"} += $percentage;
            }

            elsif ($annoType =~ /intron/){
                $typeHash{"intron"} += $percentage;
            }
            
            elsif ($annoType =~ /intergenic/){
                $typeHash{"intergenic"} += $percentage;
            }
            
            elsif( $annoType =~ /snoRNA/i
                || $annoType =~ /CDBox/i
                || $annoType =~ /SNORD/i
                || $annoType =~ /scaRNA/i
                || $annoType =~ /SNORA/i
                || $annoType =~ /HAcaBox/i
                || $annoType =~ /snRNA/i
                || $annoType =~ /miRNA/i
                || $annoType =~ /scRNA/i
                || $annoType =~ /misc_RNA/i
                || $annoType =~ /mascRNA-menRNA/i
                || $annoType =~ /ribozyme/i
                || $annoType =~ /sRNA/i 
                || $annoType =~ /sncRNA/i
                || $annoType =~ /Pol3_7SL_srpRNA/i){
                $typeHash{"sncRNA"} += $percentage;
            }

            elsif ($annoType =~ /rRNA/){
                $typeHash{"rRNA"} += $percentage;
            }           

            elsif ($annoType =~ /tRNA/){
                $typeHash{"tRNA"} += $percentage;
            }      

            elsif ($annoType =~ /lncRNA/){
                $typeHash{"lncRNA"} += $percentage;
            }      

        }

    my $outdir = "./Results/bedAnnotatorResult";
    if (! (-e $outdir)){
        executeCommand("mkdir $outdir");
    }

    my $outputFile = $outdir . "/" . $prefix . ".bedAnnotator.pie.txt";

    open ( RINPUT, ">$outputFile");
    print RINPUT "lable\t", "percent\n";

    foreach my $geneType (sort keys %typeHash) {
        my $percent = sprintf "%.2f%%", $typeHash{$geneType};
        print RINPUT "$geneType" . "_" . "$percent\t", "$typeHash{$geneType}\n";
    }

    close(RINPUT);
    close(TYPE);
}


#6 特殊的bedannotator：用位点进行注释
sub annotateKnownSites {
    my ( $confHash, $inputFile ) = @_;
    my $annotatedFile = $inputFile . ".KnownSitesAnnotator";

    my $commandLine
        = $confHash->{'bedAnnotator'} . " "##/public/home/jianhua/bioTools/bin/bedAnnotator
        . $confHash->{'bedAnnotatorParameter'}##-s 1 #1 is same strand
        . " --anno "
        . $confHash->{'splicingDir'} . "/allAnnoSites.bed6"#
        . " --bed "
        . $inputFile . " > "
        . $annotatedFile;
    &executeCommand($commandLine);

    my $KnownSitesFile  = "$outputDir/$prefix.STAR.Aligned.out.bam.startEndSeeker.KnownSitesAnnotator.KnownSites.txt";
    my $rmAnnoSitesFile = "$outputDir/$prefix.STAR.Aligned.out.bam.startEndSeeker.KnownSitesAnnotator.rmAnnoSites.txt";

    open( INPUT, "<$annotatedFile" );
    open( OUTPUT1, ">$KnownSitesFile" );
    open( OUTPUT2, ">$rmAnnoSitesFile" );

    while ( my $line = <INPUT> ) {
        $line =~ s/\s+$//;

        if (( $line =~ /#/i )) {
            print OUTPUT1 $line, "\n";
            print OUTPUT2 $line, "\n";
            next;
        }

        my @items = split /\s+/, $line;
        my $rNum = $items[6];
        my $log10PVal = $items[12];

        if (( $rNum >= 5 && $log10PVal < -4.32 )){#
            
            if (( $line =~ /intergenic/i )) {
                print OUTPUT2 $line, "\n";
            }

            else{
                print OUTPUT1 $line, "\n";
            }
        }
    }

    close(INPUT);
    close(OUTPUT1);
    close(OUTPUT2);
    return $rmAnnoSitesFile;
}

sub annotateCutSites {
    my ( $confHash, $inputFile ) = @_;
    my $annotatedFile = $inputFile . ".CutSitesAnnotator";

    my $commandLine
        = $confHash->{'bedAnnotator'} . " "##/public/home/jianhua/bioTools/bin/bedAnnotator
        . $confHash->{'bedAnnotatorParameter'}##-s 1 #1 is same strand
        . " --anno "
        . $confHash->{'splicingDir'} . "/exonStartEnd_splice5p3p_Sites.bed6"#allAnnoSites.bed6
        . " --bed "
        . $inputFile . " > "
        . $annotatedFile;
    &executeCommand($commandLine);

    my $KnownSitesFile  = "$outputDir/$prefix.STAR.Aligned.out.bam.startEndSeeker.CutSitesAnnotator.KnownCutSites.txt";
    my $rmAnnoSitesFile = "$outputDir/$prefix.STAR.Aligned.out.bam.startEndSeeker.CutSitesAnnotator.rmAnnoCutSites.txt";

    open( INPUT, "<$annotatedFile" );
    open( OUTPUT1, ">$KnownSitesFile" );
    open( OUTPUT2, ">$rmAnnoSitesFile" );

    while ( my $line = <INPUT> ) {
        $line =~ s/\s+$//;

        if (( $line =~ /#/i )) {
            print OUTPUT1 $line, "\n";
            print OUTPUT2 $line, "\n";
            next;
        }

        my @items = split /\s+/, $line;
        my $rNum = $items[6];
        my $log10PVal = $items[12];

        if (( $rNum >= 20 && $log10PVal < -4.32 )){#
            
            if (( $line =~ /intergenic/i )) {
                print OUTPUT2 $line, "\n";
            }

            else{
                print OUTPUT1 $line, "\n";
            }
        }
    }

    close(INPUT);
    close(OUTPUT1);
    close(OUTPUT2);
    return $rmAnnoSitesFile;
}


#6
sub plotPofile_old {
    my ( $confHash, $sortedBamFile, $outputDir, $prefix ) = @_;
    my $pairTag = " ";
    my $commandLine
        = $confHash->{'ncapToWig'} . " "#/public/home/huangjh/bioTools/bin/ncapToWig（没有）
        . $confHash->{'ncapToWigParameter'} . " "##-s --ncap --rpm -P --norm                #--rpm: output rpm values in wiggle file
        . $pairTag . " " . " -p "
        . $prefix . " -o "
        . $outputDir
        . " --fa "
        . $confHash->{'genomeFile'}
        . " --fai "
        . $confHash->{'faidxFile'}
        . " --bam "
        . $sortedBamFile;
    &executeCommand($commandLine);
    return 1;
}


#7
# bigWig文件是index后的二进制WIG文件，在genome browser上查看更加快速
# wigToBigWig in.wig chrom.sizes out.bw, chromsizes文件可以从UCSC上下载，就是各个染色体的长度大小
sub wigToBigwig_nap {
    my ( $confHash, $outputDir, $prefix ) = @_;
    my $style = '.count';
    if ( $confHash->{'ncapToWigParameter'} =~ /rpm/ ) {#-s --ncap --rpm -P --norm
        $style = '.rpm';
    }
    my @wiggleFiles = (
        $style . '.minus.end.wg',
        $style . '.minus.coverage.wg',
        $style . '.minus.start.wg',
        $style . '.plus.end.wg',
        $style . '.plus.coverage.wg',
        $style . '.plus.start.wg'
    );
    for ( my $i = 0; $i < scalar(@wiggleFiles); $i++ ) {
        my $wigFile = $outputDir . "/" . $prefix . $wiggleFiles[$i];
        my $bwFile  = $wigFile;
        $bwFile =~ s/\.wg/\.bw/;
        my $commandLine
            = $confHash->{'wigToBigWig'} . " "#/public/home/jianhua/softwares/kentUcscExe/wigToBigWig
            . $wigFile . " "
            . $confHash->{'chromSizeFile'} . " " . " "#/public/home/jianhua/softwares/kentUcscExe/hg38.chrom.sizes; chromsizes文件可以从UCSC上下载，就是各个染色体的长度大小
            . $bwFile;
        executeCommand($commandLine);
        $commandLine = "rm -f " . $wigFile;#画完图就删除了wig了

        executeCommand($commandLine);
    }
}

#8 
#输入文件1：剪接文件（转录起始/结束；外显子起始/结束；剪接5P/3P）
#输入文件2：四个bw文件
sub annotateTagSite {
    my ( $confHash, $outputDir, $prefix, $distSiteOutDir ) = @_;
    my $logFile
        = $distSiteOutDir . "/" . $prefix . ".annotateTagSite.log.txt";
    my $annoDir   = $confHash->{'splicingDir_2'};
    my $suffix    = "_sites.bed";
    my @annoFiles = @{ getFiles( $annoDir, $suffix ) };

    foreach my $siteFile ( sort @annoFiles ) {
        my $typeFlag = 0;
        if ( $siteFile =~ /StartSite/ || $siteFile =~ /splice5pSite/ ) {
            $typeFlag = 1;
        }
        elsif ( $siteFile =~ /EndSite/ || $siteFile =~ /splice3pSite/ ) {
            $typeFlag = 2;
        }
        my $fileType = $siteFile;
        $fileType =~ s/$suffix//;
        my $annoSiteFile = $annoDir . "/" . $siteFile;
        my $wigDir       = $outputDir . "/" . $prefix . ".rpm";
        my $startInputFile
            = $wigDir . ".plus.start.bw" . "," . $wigDir . ".minus.start.bw";
        my $endInputFile
            = $wigDir . ".plus.end.bw" . "," . $wigDir . ".minus.end.bw";
        my @wigInputFiles = ( $startInputFile, $endInputFile );

        foreach my $inputFile (@wigInputFiles) {
            my $distType = "None";
            if ( $typeFlag == 1 ) {
                $distType = "Start";
            }
            elsif ( $typeFlag == 2 ) {
                $distType = "End";
            }
            if ( $inputFile =~ /\.start\.bw/i ) {
                $distType .= "ToStart";
            }
            elsif ( $inputFile =~ /\.end\.bw/i ) {
                $distType .= "ToEnd";
            }
            my @strandParameters = ("");
            foreach my $strandParameter (@strandParameters) {
                my $siteType = $fileType;
                if ( $strandParameter =~ /anti\-strand/ ) {
                    $siteType = $siteType . "_antiStrand";
                }

                my $annotatedFile
                    = $distSiteOutDir . "/"
                    . $prefix . "."
                    . $siteType . "."
                    . $distType
                    . ".annotateSites.txt";
                #my $commandLine = "";
                my $commandLine
                    = $confHash->{'annotateSites'} . " "
                    . $confHash->{'annotateSitesParameter'} . " "
                    . $strandParameter . " "
                    . " --type "
                    . $typeFlag
                    . " --genome "
                    . $confHash->{'chromSizeFile'}
                    . " --site "
                    . $annoSiteFile
                    . " --bwg "
                    . $inputFile . " -o "
                    . $annotatedFile
                    . " 2>>$logFile ";
                &executeCommand($commandLine);

                my $annoInfo
                    = "\"Distance to "
                    . $siteType . " "
                    . $distType
                    . " sites\"";
                my $outName
                    = $distSiteOutDir . "/"
                    . $prefix . "."
                    . $siteType . "."
                    . $distType
                    . ".annotateSites";
                $commandLine
                    = " Rscript "
                    . $confHash->{'drawAnnotateSites'} . " -f "
                    . $annotatedFile . " -a "
                    . $annoInfo . " -o "
                    . $outName
                    . " 2>>$logFile ";
                &executeCommand($commandLine);
            }    # strandParameters end
        }    # input type end
    }
    return 1;
}


sub annotateBedSite {
    my ( $confHash, $outputDir, $prefix, $distSiteOutDir ) = @_;
    my $logFile
        = $distSiteOutDir . "/" . $prefix . ".annotateTagSite.log.txt";#/public/home/huangjh/NAPseq/human/NCAP_seq_condition_data/cutAdapterBarcodeData/hg38_AN-1201213_ADR_NAP-seq-rep1_star_splice/distSiteResults
    my $annoDir   = $confHash->{'RNAmodDir'};#/public/home/jianhua/genome/human/hg38/annotations/gencodeV30/annotationSite/splicing
    my $suffix    = ".bed";#_sites.bed
    my @annoFiles = @{ &getFiles( $annoDir, $suffix ) };#获取/public/home/jianhua/genome/human/hg38/annotations/gencodeV30/annotationSite/splicing路径下的_sites.bed文件($siteFile)

    foreach my $siteFile ( sort @annoFiles ) {
        my @siteFileItem = split/\./, $siteFile;
        my $siteFilePrefix = $siteFileItem[2];
        my $annoSiteFile = $annoDir . "/" . $siteFile;#用来注释的文件，还是那个剪接位置文件: _sites.bed($siteFile)（我无）


    #输出文件名；$distType = StartToStart；StartToEnd；EndToStart；EndToEnd
    #hg38_AN-2201213_ADR_NAP-seq-rep2.human_transcriptEndSite_protein_coding_rmSncRNA.EndToStart.annotateSites.txt
    my $annotatedFile
        = $distSiteOutDir . "/"
        . $prefix . "."
        . $siteFilePrefix#位点文件名
        . "ToDiffCutSite"#SiteToStart
        . ".distance.txt";

    #my $commandLine = "";

    my $commandLine
        = $confHash->{'annotateSites'} . " "#/public/home/jianhua/bioTools/bin/annotateSites（我无）
        . $confHash->{'annotateSitesParameter'} #-s(strand-specific) -m 0.0001(minimum read for each site to draw heatmap) -w 50(widow size)
        #. " " . $strandParameter . " "
        . " --type "# site type[0,1,2], 0 is coverage, 1 is start, 2 is end, default is 0
        . "0"
        . " --genome "
        . $confHash->{'chromSizeFile'}#/public/home/jianhua/softwares/kentUcscExe/hg38.chrom.sizes（我无）
        . " --site "#输入文件1：注释位点信息
        . $annoSiteFile
        . " --bed "#同时输入（正负链两个）文件2：需要注释的文件（四个bw文件）
        . $cuttingSitesBED . " > "
        . $annotatedFile#输出文件
        . " 2>>$logFile ";
    &executeCommand($commandLine);


    #R画图drawAnnotateSites.R参数
    my $annoInfo
        = "\"Distance from "
        . $siteFilePrefix . " to "#位点文件名
        . $prefix
        . " sites\"";
    my $outName
        = $distSiteOutDir . "/"
        . $prefix . "."
        . $siteFilePrefix
        . ".annotateSites";
   
   #R画图drawAnnotateSites.R
    $commandLine
        = " Rscript "
        . $confHash->{'drawAnnotateSites'} . " -f "#/public/home/jianhua/perlBioTools/annotateSite/drawAnnotateSites.R
        . $annotatedFile . " -a "#SiteToStart.annotateSites.txt
        . $annoInfo . " -o "
        . $outName
        . " 2>>$logFile ";
    &executeCommand($commandLine);

    }
}

#9
sub annotateTagDistribution {
    my ( $confHash, $inputFile, $outputDir, $prefix ) = @_;
    my $annotatedFile = $outputDir . "/" . $prefix . ".tag.bedAnnotator.txt";#/public/home/huangjh/NAPseq/human/NCAP_seq_condition_data/cutAdapterBarcodeData/hg38_AN-1201213_ADR_NAP-seq-rep1_star_splice
    my $logFile = $outputDir . "/" . $prefix . ".tag.bedAnnotator.log.txt";
    my $commandLine
        = " "
        . $confHash->{'bedAnnotator'} . " "#/public/home/jianhua/bioTools/bin/bedAnnotator
        . ' --pair --norm '
        . " --summary "#！！！！！！！总结信息，用来画饼图
        . " --anno "
        . $confHash->{'exonAnnotation'}#/public/home/jianhua/genome/human/hg38/annotations/gencodeV30/annotationBed/hg38.genecode.v30.tRNA.snoRNA.miRNA.rmsk.chrM.RefseqNcrna.curated.exonFeatures.bed6
        . " --bam "
        . $inputFile . " > "
        . $annotatedFile
        . " 2>$logFile ";
    &executeCommand($commandLine);
}

#子程序中的子程序
sub getFiles {
    my ( $dir, $suffix ) = @_;
    my @files = ();
    opendir( DIR, "$dir" ) || die "can't open the $dir directory\n";
    while ( ( my $fileName = readdir(DIR) ) ) {
        $fileName =~ s/\s+$//;
        $fileName =~ s/^\s+//;
        if ( $fileName =~ /$suffix$/ ) {
            push @files, $fileName;
        }
    }
    close(DIR);
    return \@files;
}



#没有调用
sub runStartEndSeeker {
    my ( $sortedBamFile, $platform, $readLen, $confHash ) = @_;
    my $contigFile = $sortedBamFile . ".startEndSeeker";#/public/home/jianhua/bioTools/bin/startEndSeeker

    my $pairTag    = "";
    if ( $platform eq "pairedEnd" ) {

        #$pairTag = " -P ";
    }
    my $commandLine
        = $confHash->{'startEndSeeker'} . " "#/public/home/jianhua/bioTools/bin/startEndSeeker
        . $confHash->{'startEndSeekerParameter'}#-s -t 2(minimum tag number for each start or end pos) -m 2;(minimum fold-change between start/end and nearby sites)
        . " --fa "
        . $confHash->{'genomeFile'}
        . " --fai "
        . $confHash->{'faidxFile'}
        . " --bam "
        . $sortedBamFile . " "
        . " >$contigFile";
    &executeCommand($commandLine);
    return $contigFile;
}

sub motifFinder {
    my ( $contigFile, $createDir, $prefix, $org, $extLen ) = @_;
    my %targetHash = ();
    my %seqLenHash = ();
    my $minPeakNum = 100;
    open( CONTIGRNA, "<$contigFile" ) || die "can't open the $contigFile\n";
    while ( my $line = <CONTIGRNA> ) {
        $line =~ s/\s+$//;
        if ( $line =~ /\#/ ) {
            next;
        }
        my @items = split /\s+/, $line;
        if ( scalar(@items) < 5 || $line =~ /\#/ ) {
            warn "short line: $line\n";
            next;
        }
        my $chrom      = $items[0];
        my $chromStart = $items[1];
        my $chromEnd   = $items[2];
        my $name       = $items[3];
        my $score      = $items[4];
        my $strand     = $items[5];
        $chromStart = $chromStart - $extLen;
        $chromStart = 0 if ( $chromStart < 0 );
        $chromEnd   = $chromEnd + $extLen;
        my $seqLen = $chromEnd - $chromStart;
        my $peak
            = $chrom . "\t"
            . $chromStart . "\t"
            . $chromEnd . "\t"
            . $name . "\t"
            . $score . "\t"
            . $strand;
        $seqLenHash{$name} = $seqLen;
        $targetHash{$peak} = $score;
    }    ### while
    close(CONTIGRNA);
    my $seqNum = scalar( keys %targetHash );
    #if ( $seqNum < $minPeakNum ) {
        #warn "The number of peak is less than " . $seqNum, ". Exit!!!\n";
        #return;
    #}
    my $motifDir
        = $createDir . "/" . $prefix . "_HommerGenomeResults";
    if ( !( -e $motifDir ) ) {
        &executeCommand("mkdir -p $motifDir");
    }
    my $targetFile = $motifDir . "/" . "upDownExtend" . $extLen . "nt.bed";
    my $targetFile2 = $createDir . "/" . "upDownExtend" . $extLen . "nt.bed";
    open( TARGET, ">$targetFile" );
    open( TARGET2, ">$targetFile2" );
    foreach my $peak ( sort keys %targetHash ) {
        print TARGET $peak, "\n";
        print TARGET2 $peak, "\n";
    }
    close(TARGET);
    close(TARGET2);
    my $type = "rbsSeeker";
    &findMotifs( $motifDir, $type, $prefix, $org, $targetFile );
}

sub findMotifs {
    my ( $motifDir, $type, $prefix, $org, $targetFile ) = @_;
    my $logFile = $motifDir . "/" . $prefix . "." . $type . ".log.txt";
    my $commandLine
        = "findMotifsGenome.pl "
        . $targetFile . " "
        . $org . " "
        . $motifDir . " "
        . " -norevopp -noknown -rna -len 4,5,6,7,8 -p 20 -size given -dumpFasta 2>$logFile";
    executeCommand($commandLine);
    &findMotifSites( $motifDir, $prefix, $type );
}

sub findMotifSites {
    my ( $motifDir, $prefix, $type ) = @_;
    my $tgFile     = $motifDir . "/" . "target.fa";
    my $bgFile     = $motifDir . "/" . "background.fa";
    my $motifNum   = 35;
    my $similarNum = 20;

    my $sigMotifFile
        = $motifDir . "/"
        . $prefix . "."
        . $type
        . ".significant.positionSpecifc.motifs.txt";
    my $sigMotifSeqFile
        = $motifDir . "/"
        . $prefix . "."
        . $type
        . ".significant.positionSpecifc.motifs.sequences.txt";
    my $sigMotifDistFile
        = $motifDir . "/"
        . $prefix . "."
        . $type
        . ".significant.positionSpecifc.motifs.distances.txt";
    my $logFile
        = $motifDir . "/"
        . $prefix . "."
        . $type
        . ".significant.positionSpecifc.motifs.logs.txt";

    open( SIGMOTIF, ">$sigMotifFile" )
        || die "can't open the $sigMotifFile\n";
    open( SEQMOTIF, ">$sigMotifSeqFile" )
        || die "can't open the $sigMotifSeqFile\n";
    open( DISTMOTIF, ">$sigMotifDistFile" )
        || die "can't open the $sigMotifDistFile\n";
    print SIGMOTIF
        "motifName\tmotifSeq\tdistType\tdist\tnegLg(pval)\tPercentage\tfoldChange\tdistTargetMotifNum\ttargetMotifNum\tdistBackgroundMotifNum\tbackgroundMotifNum\n";
    print SEQMOTIF
        "motifName\tmotifSeq\tdistType\tdist\tmotifInSeq\tscore\tsequence\n";
    print DISTMOTIF
        "motifName\tmotifSeq\tdistType\tdist\tnumber\tpercentage\n";

    for ( my $i = 1; $i < $motifNum; $i++ ) {
        my $motifName = "motif" . $i . ".motif";
        my $motifFile = $motifDir . "/homerResults/" . $motifName;
        my $outfile
            = $motifDir . "/" . $type . "." . $motifName . ".significant.txt";
        my $seqfile
            = $motifDir . "/"
            . $type . "."
            . $motifName
            . ".significant.seq.txt";
        my $distfile
            = $motifDir . "/"
            . $type . "."
            . $motifName
            . ".significant.dist.txt";

        if ( -e $motifFile ) {
            my $commandLine
                = "findMotif -p -10 -f 2 -s 0.85 -r 0.20 -n 20 -t $tgFile -b $bgFile -m $motifFile -o $outfile -e $seqfile -d $distfile 2>$logFile";
            executeCommand($commandLine);

            open( MOTIF, "<$outfile" ) || die "can't open the $outfile\n";
            while ( my $line = <MOTIF> ) {
                $line =~ s/\s+$//;
                if ( $line =~ /\#/ ) {
                    next;
                }
                else {
                    print SIGMOTIF $motifName, "\t", $line, "\n";
                }
            }
            close(MOTIF);

            open( SEQ, "<$seqfile" ) || die "can't open the $seqfile\n";
            while ( my $line = <SEQ> ) {
                $line =~ s/\s+$//;
                if ( $line =~ /\#/ ) {
                    next;
                }
                else {
                    print SEQMOTIF $motifName, "\t", $line, "\n";
                }
            }
            close(SEQ);

            open( DIST, "<$distfile" ) || die "can't open the $distfile\n";
            while ( my $line = <DIST> ) {
                $line =~ s/\s+$//;
                if ( $line =~ /\#/ ) {
                    next;
                }
                else {
                    print DISTMOTIF $motifName, "\t", $line, "\n";
                }
            }
            close(DIST);

            $commandLine = "rm -f $outfile $seqfile $distfile";
            executeCommand($commandLine);
        }
        for ( my $j = 1; $j < $similarNum; $j++ ) {
            $motifName = "motif" . $i . ".similar" . $j . ".motif";
            $motifFile = $motifDir . "/homerResults/" . $motifName;
            $outfile
                = $motifDir . "/"
                . $type . "."
                . $motifName
                . ".significant.txt";
            $seqfile
                = $motifDir . "/"
                . $type . "."
                . $motifName
                . ".significant.seq.txt";
            $distfile
                = $motifDir . "/"
                . $type . "."
                . $motifName
                . ".significant.dist.txt";

            if ( -e $motifFile ) {
                my $commandLine
                    = "findMotif -p -10 -f 2 -s 0.85 -r 0.20 -n 20 -t $tgFile -b $bgFile -m $motifFile -o $outfile -e $seqfile -d $distfile 2>$logFile";
                executeCommand($commandLine);
                open( MOTIF, "<$outfile" ) || die "can't open the $outfile\n";
                while ( my $line = <MOTIF> ) {
                    $line =~ s/\s+$//;
                    if ( $line =~ /\#/ ) {
                        next;
                    }
                    else {
                        print SIGMOTIF $motifName, "\t", $line, "\n";
                    }
                }
                close(MOTIF);

                open( SEQ, "<$seqfile" ) || die "can't open the $seqfile\n";
                while ( my $line = <SEQ> ) {
                    $line =~ s/\s+$//;
                    if ( $line =~ /\#/ ) {
                        next;
                    }
                    else {
                        print SEQMOTIF $motifName, "\t", $line, "\n";
                    }
                }
                close(SEQ);

                open( DIST, "<$distfile" )
                    || die "can't open the $distfile\n";
                while ( my $line = <DIST> ) {
                    $line =~ s/\s+$//;
                    if ( $line =~ /\#/ ) {
                        next;
                    }
                    else {
                        print DISTMOTIF $motifName, "\t", $line, "\n";
                    }
                }
                close(DIST);

                $commandLine = "rm -f $outfile $seqfile $distfile";
                executeCommand($commandLine);
            }
        }
    }
    close(SIGMOTIF);
    close(SEQMOTIF);
    close(DISTMOTIF);
}