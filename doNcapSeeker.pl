#!/usr/bin/perl -w

BEGIN {
    unshift( @INC, "/public/home/huangjh/perlBioTools/packages" );
}

use Getopt::Long;
use strict;
use commonTools;

my $sampFile = "";
my $sufxName = ".clip5p3pAdapter.pair.fq.gz";
my $progDir  = "/public/home/huangjh/perlBioTools/ncapSeeker";
my $help     = "";
if ( scalar(@ARGV) < 1 ) {
    usage();
}

my $optLong = GetOptions(
    "samp=s" => \$sampFile,    # configure file
    "sufx=s" => \$sufxName,    # suffix Name of file
    "prog=s" => \$progDir,     # program Directory
    "help|h" => \$help
) or usage();

&usage() if $help;

sub usage
### usage function
{
    print
        "Usage: perl $0 [options] --samp <sample file> --sufx <suffix name>\n";
    print "--samp        sample file\n";
    print "--sufx        suffix name of files[default=.clip.pair.fq.gz]\n";#实际是.clip5p3pAdapter.pair.fq.gz
    print
        "--prog        program Directory[default=/public/home/jianhua/perlBioTools/ncapSeeker]\n";
    print "--help        help information\n";
    exit;
}

$sufxName = ".cutadapt.collapse.cutBarcodes.fq.gz" if ($sampFile~~/human_RIP_NAP_seq_merge_samples/);
&runAllRNAdata( $sampFile, $sufxName, $progDir );

sub runAllRNAdata {
    my ( $infoFile, $sufxName, $progDir ) = @_;
    my $i = 1;
    open( INFO, "<$infoFile" ) || die "can't open the $infoFile\n";
    while ( my $line = <INFO> ) {
        $line =~ s/\s+$//;
        next if ( $line =~ /genome/ );
        my @items     = split /\t/, $line;
        my $org       = $items[0];
        my $srp       = $items[1];
        my $ftp       = $items[2];
        my $sample    = $items[3];
        my $barcode   = $items[4];
        my $readLen   = $items[5];
        my $seqType   = $items[6];
        my $sampleType   = $items[7];
        my $platform  = "pairedEnd";
        my $seqFormat = "fastq";
        my @sraInfo   = split /\//, $ftp;
        my $sra       = $sraInfo[$#sraInfo];
        my $createDir = "./" . $srp . "/" . $org . "_" . $sra . "_" . $sample . "_star_splice";


        my $sra1File  = "./" . "$srp" . "/" . $sra . "_1" . $sufxName;
        my $sra2File  = "./" . "$srp" . "/" . $sra . "_2" . $sufxName;

        if ( $sample ~~ /SHAPE/i ) {
            $createDir = "./" . $srp . "/" . $org . "_" . $sra . "_" . $sample . "_star_splice";
            $sufxName=".clip5p3pAdapter.pair.fastq.gz";
            $sra1File  = "./" . "$srp" . "/" . $sra . "_R1" . $sufxName;
            $sra2File  = "./" . "$srp" . "/" . $sra . "_R2" . $sufxName;
        }

        if ( $sample ~~ /DMSO|NAI-N3/i ) {
            $createDir = "./" . $srp . "/" . $org . "_" . $sra . "_" . $sample . "_star_alignToGenome";
            $sufxName=".fastq.gz";
            $sra1File  = "./" . "$srp" . "/" . $sra . "_R1" . $sufxName;
            $sra2File  = "./" . "$srp" . "/" . $sra . "_R2" . $sufxName;
        }


        if ( $sample ~~ /nanopore/i ) {
            $createDir = "./" . $srp . "/" . $org . "_" . $sra . "_" . $sample . "_minimap_splice";
        }

        if ( $sample ~~ /RNAseq/ ) {
            $sufxName=".cutadapt.fq.gz";
            $platform  = "pairedEnd";
            $seqFormat = "fastq";
            $sra1File  = "./" . "cutAdapterBarcodeData" . "/" . $sra . "_R1" . $sufxName;
            $sra2File  = "./" . "cutAdapterBarcodeData" . "/" . $sra . "_R2" . $sufxName;
        }


        if ( $sampleType eq "treat" ) {
            $platform  = "pairedEnd";
            $seqFormat = "fastq";
            $sra1File  = "./" . "NCAP_seq_condition_data/cutAdapterBarcodeData" . "/" . $sra . "_1" . $sufxName;
            $sra2File  = "./" . "NCAP_seq_condition_data/cutAdapterBarcodeData" . "/" . $sra . "_2" . $sufxName;
        }

        elsif ( $sampleType eq "normal" ) {
            $platform  = "pairedEnd";
            $seqFormat = "fastq";
            $sra1File  = "./" . "cutAdapterBarcodeData" . "/" . $sra . "_1" . $sufxName;
            $sra2File  = "./" . "cutAdapterBarcodeData" . "/" . $sra . "_2" . $sufxName;
        }
        if ( !( -e $createDir ) ) {
            &executeCommand("mkdir $createDir");#cutAdapterBarcodeData
        }
        my $prefix  = $org . "_" . $sra . "_" . $sample;
        my $logFile = $createDir . "/" . $sra . ".all.ncapSeeker.log.txt";
        my $commandLine
            = "perl "
            . $progDir
            . "/ncapSeeker.pl --conf "
            . $progDir
            . "/ncapSeekerConfigure_"
            . $org . ".xml "
            . " --output "
            . $createDir
            . " --prefix "
            . $prefix
            . " --platform "
            . $platform
            . " --format "
            . $seqFormat
            . " --readLen "
            . $readLen
            . " --read1 "
            . $sra1File
            . " --org "
            . $org
            . " 2>>$logFile &";
        if ( $platform eq "pairedEnd" ) {
            $commandLine
                = "perl "
                . $progDir
                . "/ncapSeeker.pl --conf "
                . $progDir
                . "/ncapSeekerConfigure_"
                . $org . ".xml "
                . " --output "
                . $createDir
                . " --prefix "
                . $prefix
                . " --platform "
                . $platform
                . " --format "
                . $seqFormat
                . " --readLen "
                . $readLen
                . " --read1 "
                . $sra1File
                . " --read2 "
                . $sra2File
                . " --org "
                . $org
                . " 2>>$logFile &";
        }
        executeCommand($commandLine);
        &limitProgramNum( "ncapSeeker.pl", 20 );
        $i++;
    }
    close(INFO);
}
