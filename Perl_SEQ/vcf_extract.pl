#!/usr/bin/perl

use strict;
use warnings;
use Getopt::Long;
use List::Util qw (first);

# Initialize the parameters
my ($help, $by, $file, $list);
my $info = "ALL";
# Read the parameters
usage() if ( @ARGV < 1 or 
	     ! GetOptions('help|?' => \$help,
			  'by=s'   => \$by,
			  'file=s' => \$file,
			  'list=s' => \$list,
			  'info=s' => \$info) or
	     defined $help );

# Configure list (individual or snps) 
my $fh1_content = do {
    local $/;
    local @ARGV = ( $list );
    <>
};


my (@indList, @indIndex);
my (@snpList, %chr_pos);

if ($by =~ /ind/i) {
    @indList = split "\n", $fh1_content;
} elsif ($by =~ /snp/i) {
    @snpList = split "\n", $fh1_content;
    %chr_pos = map {$_ => 1} @snpList;
}

# Process vcf files
open my $fh, "<$file" or die "Cannot open vcf file: $file";
my @header;
my %hash;

while ( <$fh> ) {
    chomp; 
    if ( /^##\s?\w+/ ) { # header beginning with two 

	next;

    } elsif ( /^#\s?\w+/ ) { # header containing individuals

	@header = split "\t", $_;

	if ( $by =~ /ind/ ) {
	    foreach my $sample (@indList) {
		my $id = first { $header[$_] eq $sample } 9..$#header;
		push (@indIndex, $id);
	    } # end of foreach
	} else {
	    @indIndex = (9..$#header);
	} # end of if-else
	 
    }  else {

	my @row = split /\t/, $_;
	my $chr_pos = $row[0] . "_" . $row[1];
	my @selected_ind;
	if  ( $by =~ /ind/i or ($by =~ /snp/i && exists $chr_pos{$chr_pos}) ) { 
	    print join ("\t",@row[0..8]) , "\t";
	    my @format = split ":", $row[8];
	    foreach my $i (@indIndex) {
		if ( $row[$i] !~ /^\.$/ ) {
		    @{ $hash{$header[$i]} }{@format} = split ":", $row[$i]; # hash slice @{ $hash{$key} } {@keys}
		    ($info=~/ALL/i) ? push (@selected_ind, $row[$i]) : push (@selected_ind, $hash{$header[$i]}{$info});
		} elsif ( $row[$i] eq "." ) {
		    my $length = scalar @format;
		    my @uncalled = map { "\." } 1..$length;
		    @{ $hash{$header[$i]} }{@format}  = @uncalled;
		    ($info=~/ALL/i) ? push (@selected_ind, $row[$i]) : push (@selected_ind, $hash{$header[$i]}{$info});
		} # end of if-els-if
		
	    }# end of foreach
	    print join("\t",@selected_ind), "\n";
	} else {
	    next;
	} # end of if-else

    } # end of if-elsif-else:header-header-snp
} # end of while
	    
sub usage {
    print "Unkonwn options: @_\n" if (@_);
    print "\nvcf_extrac.pl\n";
    print "2015\n";
    print "By Xiaodong Liu \t\t xiaodong.liu\@ebc.uu.se\n";
    print "###############################################\n";
    print "usage: vcf_extract.pl [--by ind(individuals) or snp(SNPs)] [--file vcf][--list file] [--help or -?]\n";
    print "###############################################\n\n";
    exit ;
}
