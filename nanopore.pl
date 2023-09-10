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
use XML::Simple;
use File::Basename;

$| = 1;

my $configureXmlFile = "xml";
my $readFile         = "";
my $outputDir        = "./output";
my $prefix           = "prefix";
my $platform         = "singleEnd";
my $readLen          = 1500;
my $seqFormat        = "fastq";
my $help             = 0;
my $org              = "hg38";

if ( scalar(@ARGV) < 4 ) {
    usage();
}

my $optLong = GetOptions(
    "conf=s"     => \$configureXmlFile,    # configure file
    "read=s"     => \$readFile,            # read1File file
    "output=s"   => \$outputDir,           # outputDir
    "prefix=s"   => \$prefix,              # prefix
    "platform=s" => \$platform,            # platform
    "readLen=s"  => \$readLen,             # readLen
    "format=s"   => \$seqFormat,           # sequence format
    "org=s"      => \$org,                 # org
    "help|h"     => \$help
) or usage();

&usage() if $help;

### configure files
my $xmlConfHash = XMLin($configureXmlFile);

my $shell="
cd /public/home/huangjh/NAPseq/nanopore_3NGS && \
nohup perl /public/home/huangjh/perlBioTools/ncapSeeker/codes/nanopore/doNanopore.pl \
--samp mouse_NcapNanopore_C2C12_samples.txt &
";

sub usage
### usage function
{
    print
        "Usage: perl $0 [options] --conf <configure file> --output <output file> --prefix <file prefix> --read <read file> --org <organism, e.g. hg38> \n";
    print "--conf        conf file\n";
    print "--read       read  file\n";
    print "--output      output dir name[default=output]\n";
    print "--prefix      prefix for file\n";
    print
        "--platform    platform [default=singleEnd], pairedEnd or singleEnd\n";
    print "--readLen     read length [default=101]\n";
    print "--format      seq format [default=fastq], fastq or fasta\n";
    print "--org         organism [default=hg38]\n";
    print "--help        help information\n";
    exit;
}

if ( !( -e $outputDir ) ) {
    warn "ceate output dir\n";
    &executeCommand("mkdir $outputDir");
}

my $sortedBamFile = &alignNanoporeReadToGenome( $xmlConfHash, $readFile, $outputDir, $prefix );
my $sortedBamFile = $outputDir . "/" . $prefix . ".minimap2.sorted.bam";

##################################################################################################################
my $bamFile = $sortedBamFile;
### 1.cal how many read containing p5-p3-p8 adapters.
# my $readSummaryFile = &computeStartEndTag($xmlConfHash, $readFile, $outputDir, $prefix);

## 2.correlation of known gene expression cal by bamExpression.
# my $bamExpression = &runBamExpression( $xmlConfHash, $xmlConfHash->{'refSeqBed'}, $bamFile, $outputDir, $prefix );

## 3.library tag distribution.
# my $annotateTagInfoFile = $outputDir . "/" . $prefix . ".tag.bedAnnotator.txt";
# if (!(-e $annotateTagInfoFile)){
#     my $annotateTagInfoFile = &annotateTagDistribution( $xmlConfHash, $bamFile, $outputDir, $prefix );
# }
# my $prefix_short="";
# my $seqType="CAN-seq";
# $seqType="RNA-seq" if ($prefix=~/RNAseq/i);
# if ($prefix=~/diff(\d+)d.*rep(\d+)/i){
#     $prefix_short="TGS-$seqType-diff$1d-rep$2";
# }
# &pieplotTagInfo( $xmlConfHash, $annotateTagInfoFile, $outputDir, 0, 14, $prefix_short);
##################################################################################################################


#### fetch nanopore full-length RNA:
# my $bamtobed = &bamtobed($sortedBamFile);
# my $bamtobed = $outputDir . "/" . $prefix . ".minimap2.sorted.bam.bamtobed";
# my $bamtobed_uniq = &uniq( $bamtobed );

# my $canSeqFile = &runCanSeeker( $xmlConfHash, $sortedBamFile, $outputDir, $prefix );
# my $evoConGenomeFile     = &getResultConVals( $xmlConfHash, $canSeqFile);
# my $annoteGenomeResultFile  = &annotateResults( $xmlConfHash, $evoConGenomeFile);

# my $plot = &plotPofile( $xmlConfHash, $sortedBamFile, $outputDir, $prefix );
# my $wiggle = &wigToBigwig( $xmlConfHash, $outputDir, $prefix );


# my $distSiteOutDir = $outputDir . "/" . "distSiteResults";
# if ( !( -e $distSiteOutDir ) ) {
#     warn "ceate distSiteOutDir dir\n";
#     &executeCommand("mkdir $distSiteOutDir");
# }
# &annotateTagSite( $xmlConfHash, $outputDir, $prefix, $distSiteOutDir );

#&annotateTagDistribution( $xmlConfHash, $bamFile, $outputDir, $prefix );

sub runBamExpression {
    my ( $confHash, $annoFile, $bamFile, $outputDir, $prefix ) = @_;
    my $outfile   = $outputDir . "/" . $prefix . ".bamExpression.txt";
    my $logFile   = $outputDir . "/" . $prefix . ".bamExpression.log.txt";
    my $commandLine
        = $confHash->{'bamExpression'} . " "
        . $confHash->{'bamExpressionParameter'} . " "#--pair(paired-end) --noSplice(discard reads span exon-exon) --norm
        . "--anno " . $annoFile
        . " --bam " . $bamFile
        . " -o "    . $outfile
        . " 2>$logFile";
    executeCommand($commandLine);
    return $outfile;
}

sub pieplotTagInfo {
    my ( $confHash, $summaryFile, $outDir, $annoIdx, $plotTopNums,$prefix ) = @_;

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
        my $percet = $items[2];

        if( $annoType=~/Pol3_7SL_srpRNA/i ){
            $typeHash{"7SLRNA"}+=$percet;
            $sum+=$percet; 
        }

        elsif( $annoType=~/rRNA/i ){
            $typeHash{"rRNA"}+=$percet;
            $sum+=$percet; 
        }

        elsif( $annoType=~/snRNA/i ){
            $typeHash{"snRNA"}+=$percet;
            $sum+=$percet; 
        }

        elsif( $annoType=~/snoRNA/i || $annoType=~/snoRd/i || $annoType=~/snoRa/i ){
            $typeHash{"snoRNA"}+=$percet;
            $sum+=$percet;
        }

        elsif( $annoType=~/piRNA/i ){
            $typeHash{"piRNA"}+=$percet;
            $sum+=$percet;
        }

        elsif( $annoType=~/miRNA/i ){
            $typeHash{"miRNA"}+=$percet;
            $sum+=$percet;
        }

        elsif( $annoType=~/cds/i){
            $typeHash{"mRNA"}+=$percet;
            $sum+=$percet;
        }

        elsif( $annoType=~/utr3/i ){
            $typeHash{"mRNA"}+=$percet;
            $sum+=$percet;
        }

        elsif( $annoType=~/utr5/i ){
            $typeHash{"mRNA"}+=$percet;
            $sum+=$percet;
        }

        elsif( $annoType=~/intron/i ){
            $typeHash{"mRNA"}+=$percet;
            $sum+=$percet;
        }   

        elsif( $annoType=~/lncRNA/i ){
            $typeHash{"lncRNA"}+=$percet;
            $sum+=$percet;
        } 

        elsif( $annoType=~/snRNA/i ){
            $typeHash{"snRNA"}+=$percet;
            $sum+=$percet;
        }        

        elsif( $annoType=~/repeat/i ){
            $typeHash{"repetitive-elements"}+=$percet;
            $sum+=$percet;
        }

        elsif( $annoType=~/intergenic/i ){
            $typeHash{"intergenic"}+=$percet;
            $sum+=$percet;
        }
        elsif( $annoType=~/pseudogene/i ){
            $typeHash{"pseudogene"}+=$percet;
            $sum+=$percet;
        }
        else{
            $typeHash{"others"}+=$percet;
            $sum+=$percet;
        }
    }

    my $calculator =0;
    my @plot_arry=();
    foreach my $key (sort { $typeHash{$b}<=> $typeHash{$a} } keys %typeHash){
        if ($calculator<$plotTopNums){
            # my $percent = sprintf "%.10f", $typeHash{$key}/$sum;
            my $percent = $typeHash{$key};
            $percent = sprintf("%.2f", $percent);
            print OUT "$key\t", "$percent\n";
            push @plot_arry,$key;
        }
        $calculator++;
    }
    
    close(TYPE);
    close(OUT);
    undef %typeHash;
    executeCommand("Rscript /public/home/huangjh/perlBioTools/canSeq/pieplot_tag.R -f $outFile -o $outFile -p $prefix");

}


sub computeStartEndTag {
    my ( $confHash, $readFile, $outputDir, $prefix ) = @_;
    my $outFile = $outputDir . "/" . $prefix . ".startEndTagSummary.txt";
    open( OUTPUT, ">$outFile" ) || die "can't open the $outFile\n";
    print OUTPUT "#sampleName", "\t", "totalRead", "\t", "p3Read", "\t", "p5Read", "\t", "p8Read", "\t", "othRead", "\t", "startEnds", "\t","percentage", "\n";
    warn "# deal with the file: " . $readFile, "\n";
    my $totalRead = 0;
    my $p3Read    = 0;
    my $p5Read    = 0;
    my $p8Read    = 0;
    my $othRead   = 0;
    open( my $rfh, "<$readFile" ) || die("Could not open gzipped read file: $readFile \n");

    while ( my $line = <$rfh> ) {
        $line =~ s/\s+$//;
        if ($line~~/runid/){
           $totalRead++;
            if ( $line =~ /\:3\s+runid/ ) {
                $p3Read++;
            }
            elsif ( $line =~ /\:5\s+runid/ ) {
                $p5Read++;
            }
            elsif ( $line =~ /\:8\s+runid/ ) {
                $p8Read++;
            }
            else {
                $othRead++;
            }         
        }
    }
    close($rfh);
    print OUTPUT $readFile, "\t", $totalRead, "\t", $p3Read, "\t",$p5Read, "\t",$p8Read, "\t", $othRead, "\t", ( $p3Read + $p5Read + $p8Read ), "\t", ( $p3Read + $p5Read + $p8Read ) / $totalRead * 100,"\n";
    close(OUTPUT);
    return $outFile;
}


sub bamtobed{
    my ($bam) = @_;
    my $bamtobed = $bam . ".bamtobed";
    executeCommand("bedtools bamtobed -i $bam -bed12 > $bamtobed &");
    return $bamtobed;
}

sub uniq {
    my ( $inputFile ) = @_;
    open( IN, "<$inputFile" ) || die "can't open the $inputFile\n";
    my $outFile = $inputFile;
    $outFile =~ s/\.txt$//;
    $outFile =~ s/\.bed$//;
    $outFile .= ".with358Adapters.uniq.bed";
    my $bed12_header = "#chrom\tchromStart\tchromEnd\tname\tlength\tstrand\tthickStart\tthickEnd\titemRgb\tblockCount\tblockSizes\tblockStarts";
    open( OUT, ">$outFile" ) || die "can't open the $outFile\n";
    print OUT "$bed12_header\tcounts\n";    
    my %lineHash=();
    my $nameHash=();
    my %countHash=();
    while ( my $line = <IN> ) {
        $line =~ s/\s+$//;
        if ($line ~~/#chrom/i){
            next;
        }
        my @items  = split /\t/, $line;
        my $chrom      = $items[0];
        my $chromStart = $items[1];
        my $chromEnd   = $items[2];
        my $name       = $items[3];
        my $score      = $items[4];
        my $strand     = $items[5];
        my $extendSeq  = $items[6];
        next if ($name!~/:[358]/);
        my $RNA = $chrom."_".$chromStart."_".$chromEnd."_".$strand;### input as $a, $b
        $countHash{$RNA}++;
        $lineHash{$RNA}=$line;
    }

    foreach my $key ( sort by_chrom keys %countHash){
        my $line=$lineHash{$key};
        my @items  = split /\t/, $line;
        my $lineInfojoin1  = join "\t", @items[0..3];
        my $length  = $items[2]-$items[1];
        my $lineInfojoin2  = join "\t", @items[5..$#items];
        print OUT $lineInfojoin1,"\t",$length,"\t",$lineInfojoin2,"\t",$countHash{$key},"\n";
    }
    close(IN);
    close(OUT);
    return $outFile;
}

sub by_chrom {
    my $chr_a=(split/_/,$a)[0];
    $chr_a=~s/chr//i;
    my $chrStart_a=(split/_/,$a)[1];
    my $chrEnd_a=(split/_/,$a)[2];

    my $chr_b=(split/_/,$b)[0];
    $chr_b=~s/chr//i;
    my $chrStart_b=(split/_/,$b)[1];
    my $chrEnd_b=(split/_/,$b)[2];

    $chr_a cmp $chr_b or $chrStart_a<=>$chrStart_b or $chrEnd_a<=>$chrEnd_b;
}

sub getResultConVals {
    my ( $confHash, $inputFile ) = @_;
    my $conValFile = &getEvoConVals(
        $confHash->{'bedConservation'},
        $confHash->{'bedConservationParameter'},
        $confHash->{'conservationProfile'}, $inputFile
    );
    return $conValFile;
}

sub annotateResults {
    my ( $confHash, $inputFile ) = @_;
    my $annotatedFile = $inputFile . ".bedAnnotator.txt";
    my $commandLine
        = $confHash->{'bedAnnotator'} . " "
        . $confHash->{'bedAnnotatorParameter'}
        . " --anno "
        . $confHash->{'exonAnnotation'}
        . " --bed "
        . $inputFile . " > "
        . $annotatedFile;
    &executeCommand($commandLine);
    return $annotatedFile;
}

sub annotateTagDistribution {
    my ( $confHash, $inputFile, $outputDir, $prefix ) = @_;
    my $annotatedFile = $outputDir . "/" . $prefix . ".tag.bedAnnotator.txt";
    my $logFile = $outputDir . "/" . $prefix . ".tag.bedAnnotator.log.txt";
    my $commandLine
        = " "
        . $confHash->{'bedAnnotator'} . " "
        . '  -s 1 --norm '
        . " --summary "
        . " --anno "
        . $confHash->{'exonAnnotation'}
        . " --bam "
        . $inputFile . " > "
        . $annotatedFile
        . " 2>$logFile ";
    &executeCommand($commandLine);
}

sub annotateTagSite {
    my ( $confHash, $outputDir, $prefix, $distSiteOutDir ) = @_;
    my $logFile = $distSiteOutDir . "/" . $prefix . ".annotateTagSite.log.txt";
    # my $annoDir   = $confHash->{'splicingDir'};
    # my $annoDir   = $confHash->{'commonDir'};
    my $annoDir   = $confHash->{'genecodeDir'};
    my $suffix    = "_sites.bed";
    # my $suffix    = "1nt_sites.bed";
    my @annoFiles = @{ getFiles( $annoDir, $suffix ) };

    foreach my $siteFile ( sort @annoFiles ) {
        next if ($siteFile!~/HIST/);
        my $scaleFactor=10;
        if ($siteFile~~/HIST/){
            $scaleFactor=0.01;
        }
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
                    . $inputFile . " > "
                    . $annotatedFile
                    . " 2>>$logFile ";
                # &executeCommand($commandLine);

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
                    . $confHash->{'drawAnnotateSites'} 
                    . " -f " . $annotatedFile 
                    . " -a " . $annoInfo 
                    . " -o ". $outName
                    . " -s ". $scaleFactor
                    . " 2>>$logFile ";
                &executeCommand("$commandLine &");
            }    # strandParameters end
        }    # input type end
    }
    return 1;
}


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


sub runCanSeeker {
    my ( $confHash, $sortedBamFile, $outputDir, $prefix ) = @_;
    my $canSeeker = $outputDir . "/" . $prefix . ".genome.canSeeker";
    my $commandLine
        = $confHash->{'canSeeker'} . " "
        . $confHash->{'canSeekerParameter'} . " "
        . " --fa "
        . $confHash->{'genomeFile'}
        . " --fai "
        . $confHash->{'faidxFile'}
        . " --input "
        . $sortedBamFile . " "
        . " >$canSeeker";
    &executeCommand($commandLine);
    return $canSeeker;
}

sub alignNanoporeReadToGenome {
    my ( $confHash, $readFile, $outputDir, $prefix ) = @_;

    my $samOutFile = $outputDir . "/" . $prefix . ".minimap2.sam";
    my $commandLine
        = $confHash->{'minimap'} . " "
        . " --junc-bed "
        . $confHash->{'junctionBed'} . " "
        . $confHash->{'minimapParameter'} . " -o "
        . $samOutFile . " "
        . $confHash->{'minimapIdx'} . " "
        . $readFile;
    &executeCommand($commandLine);

    my $bamOutFile = $outputDir . "/" . $prefix . ".minimap2.bam";
    &samToBam( $confHash->{'samtools'}, $samOutFile, $bamOutFile );
    my $sortBamFile = $outputDir . "/" . $prefix . ".minimap2.sorted.bam";
    &sortBam( $confHash->{'samtools'}, $bamOutFile, $sortBamFile );
    &indexBam( $confHash->{'samtools'}, $sortBamFile );

    $commandLine = " rm -f $bamOutFile $samOutFile";
    &executeCommand($commandLine);

    return ($sortBamFile);
}

sub plotPofile {
    my ( $confHash, $sortedBamFile, $outputDir, $prefix ) = @_;
    my $pairTag = " ";
    my $commandLine
        = $confHash->{'ncapToWig'} . " "
        . $confHash->{'ncapToWigParameter'} . " "
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

sub wigToBigwig {
    my ( $confHash, $outputDir, $prefix ) = @_;
    my $style = '.count';
    if ( $confHash->{'ncapToWigParameter'} =~ /rpm/ ) {
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
