#!/usr/bin/perl -w

BEGIN {
    unshift( @INC, "/public/home/huangjh/perlBioTools/packages" );
}

use Getopt::Long;
use strict;
use commonTools;

my $sampFile = "";
my $sufxName = ".cutadapt.all.fq";
my $progDir  = "/public/home/huangjh/perlBioTools/ncapSeeker/codes/nanopore";
my $help     = "";
if ( scalar(@ARGV) < 2 ) {
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
    print "--sufx        suffix name of files[default=.clip.pair.fq.gz]\n";
    print
        "--prog        program Directory[default=/public/home/huangjh/perlBioTools/canSeqNanopore]\n";
    print "--help        help information\n";
    exit;
}

&runAllRNAdata( $sampFile, $sufxName, $progDir );

sub runAllRNAdata {
    my ( $infoFile, $sufxName, $progDir ) = @_;
    my $i = 1;
    open( INFO, "<$infoFile" ) || die "can't open the $infoFile\n";
    while ( my $line = <INFO> ) {
        $line =~ s/\s+$//;
        next if ($line =~ /genome/);
        my @items     = split /\t/, $line;
        my $org       = $items[0];
        my $srp       = $items[1];
        my $sra       = $items[2];
        my $sample    = $items[3];
        my $barcode   = $items[4];
        my $readLen   = $items[5];
        my $seqType   = $items[6];
        my $platform  = "singleEnd";
        my $seqFormat = "fasta";

        my $createDir
            = "./" . $srp . "/" . $org . "_" . $sra . "_" . $sample . "_minimap_splice";
        my $sraFile = "./" . $srp . "/" . $sra . $sufxName;
        if ( !( -e $sraFile ) ) {
            warn "is not exist ", $sraFile, "\n";
        }
        else {
            if ( !( -e $createDir ) ) {
                &executeCommand("mkdir $createDir");
            }
            my $prefix  = $org . "_" . $sra . "_" . $sample."_splice";
            my $logFile = $createDir . "/" . $sra . ".all.nanopore.log.txt";
            my $commandLine
                = "nohup perl "
                . $progDir
                . "/nanopore.pl --conf "
                . $progDir
                . "/nanoporeConfigure_"
                . $org . ".xml "
                . " --output "
                . $createDir
                . " --prefix "
                . $prefix
                . " --platform "
                . $platform
                . " --format "
                . $seqFormat
                . " --read "
                . $sraFile
                . " --org "
                . $org
                . " 2>>$logFile &";
            executeCommand($commandLine);
            &limitProgramNum( "nanopore.pl", 1 );
            $i++;
        }
    }
    close(INFO);
}
