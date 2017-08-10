#!/usr/bin/perl

use strict;
use warnings;
use lib './lib';
use vcf2qtl qw ( open_input uniq phase_check phase transform first );

my $infile = shift;
my $fh = open_input($infile);
my %chr_pos;
my %hash;
my ( @header, @line_f1f_index, @line_f1m_index, @f2prog_index);
my ($f1f, $f0_f1f_f, $f0_f1f_m) = ("33C1", "354P", "24_2B"); #female line (f refers to female)
my ($f1m, $f0_f1m_f, $f0_f1m_m) = ("18R1", "205B", "303P"); # male line (m refers to male)


while (<$fh>) {
    chomp;
    if ( /^##\s?\w+/ ) { # header beginning with two #
	next;

    } elsif ( /^#\s?\w+/ ) { # header contaning names
	@header = split "\t", $_;
	
	foreach my $sample ($f1f, $f0_f1f_f,$f0_f1f_m) {
	    my $id = first { $header[$_] eq $sample } 9..$#header;
	    die "No match $sample" if ( length($id) <=0 );
	    push @line_f1f_index, $id; # need to be stored in order, otherwise progeny, female and male will be mixed
	}

	foreach my $sample ($f1m, $f0_f1m_f, $f0_f1m_m) {
	    my $id = first { $header[$_] eq $sample } 9..$#header;
	    die "No match $sample" if ( length($id) <=0 );
	    push @line_f1m_index, $id; # need to be stored in order, otherwise progeny, female and male will be mixed 
	}
	
	@f2prog_index = grep {$header[$_] =~ /^G[0-9]+/} 9..$#header;
	print "Marker\t", join("\t",@header[@f2prog_index]), "\n";

    } else { #### snps
	my @row = split /\t/, $_;
	my @format = split ":", $row[8];
	my (@line_f1f_GT, @line_f1m_GT, @f2prog_GT);
	my (@phased_f1f, @phased_f1m, @phased_f2prog);

	{ ### phase line f1 female
	print $row[0] . "_" . "$row[1]","\t";    
	foreach my $i (@line_f1f_index) {
	    
	    if ( $row[$i] !~ /^\.$/ ) {
		@{ $hash{$header[$i]} }{@format} = split ":", $row[$i]; # hash slice @{ $hash{$key} } {@keys}
		print "$hash{$header[$i]}{'GT'}\t";
	    } elsif ( $row[$i] eq "." ) {
		my $length = scalar @format;
		my @uncalled = map {"\."} 1..$length;
		@{ $hash{$header[$i]} }{@format} = @uncalled;
		print "$hash{$header[$i]}{'GT'}\t";
		
	    }
	    push @line_f1f_GT, $hash{$header[$i]} {'GT'};

	} # end of foreach i
    
	my @prog_line_f1f = @line_f1f_GT[0..($#line_f1f_GT - 2)]; # f1f
	my $prog_ref_line_f1f = \@prog_line_f1f;
	@phased_f1f = phase ( $prog_ref_line_f1f, $line_f1f_GT[$#line_f1f_GT - 1], $line_f1f_GT[$#line_f1f_GT] );
	print "@phased_f1f\n";
	} ### end of phase line f1 female
	
	{ ### phase line f1 male
	print $row[0] . "_" . "$row[1]","\t";
	foreach my $j (@line_f1m_index) {
	    
	    if ( $row[$j] !~ /^\.$/ ) {
		@{ $hash{$header[$j]} }{@format} = split ":", $row[$j]; # hash slice @{ $hash{$key} } {@keys}
		print "$hash{$header[$j]}{'GT'}\t";
	    } elsif ( $row[$j] eq "." ) {
		my $length = scalar @format;
		my @uncalled = map {"\."} 1..$length;
		@{ $hash{$header[$j]} }{@format} = @uncalled;
		print "$hash{$header[$j]}{'GT'}\t";
		
	    }
	    push @line_f1m_GT, $hash{$header[$j]} {'GT'};

	} # end of foreach j

	my @prog_line_f1m = @line_f1m_GT[0..($#line_f1m_GT - 2)]; # f1m
	my $prog_ref_line_f1m = \@prog_line_f1m;
	@phased_f1m = phase ( $prog_ref_line_f1m, $line_f1m_GT[$#line_f1m_GT - 1], $line_f1m_GT[$#line_f1m_GT] );
	print "@phased_f1m\n";
	} ### end of phase line f1 male

	{ ### phase F2 progenies
	foreach my $k (@f2prog_index) {
	    
	    if ( $row[$k] !~ /^\.$/ ) {
		@{ $hash{$header[$k]} }{@format} = split ":", $row[$k]; # hash slice @{ $hash{$key} } {@keys}                                       
  		#print "$hash{$header[$k]}{'GT'}\t";
            } elsif ( $row[$k] eq "." ) {
		my $length = scalar @format;
                my @uncalled = map {"\."} 1..$length;
                @{ $hash{$header[$k]} }{@format} = @uncalled;
                #print "$hash{$header[$k]}{'GT'}\t";

            }
            push @f2prog_GT, $hash{$header[$k]} {'GT'};

	} # end of foreach k
	
	my $f2prog_GT_ref = \@f2prog_GT;
	@phased_f2prog = phase ($f2prog_GT_ref, $phased_f1f[0], $phased_f1m[0]);
	print $row[0] . "_" . "$row[1]","\t";
	print join ("\t", @phased_f2prog),"\n";
	} ### end of phase line F2 progenies

	{ ### transformation
	my $phased_f2prog_ref = \@phased_f2prog;
	my @QTLcode = transform ($phased_f2prog_ref, @phased_f1f, @phased_f1m);
	print $row[0] . "_" . "$row[1]","\t";
	print join("\t", @QTLcode), "\n";
	} ### end of transformation



    } #### end of if-elsif-else

} # end of while
