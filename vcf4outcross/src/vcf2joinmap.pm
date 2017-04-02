package vcf2joinmap;

use strict;
use warnings;

use PerlIO::gzip;
use Env qw ( @PATH );
use List::Util qw ( first max );
use List::Util qw (first);
use Data::Dumper;
use feature qw(switch);

use Exporter qw(import);
our @EXPORT_OK = qw(open_input check_GT_pedigree segCodeConvert segCodeCheck);

sub open_input {
    my $filename = shift;
    my $fh;
    if ( ! -f $filename ) {
        print "Cannot open file: $filename\n";
        exit;
    } else {
        if ($filename =~ /\.gz$/) {
            open $fh, "<:gzip", $filename;
        } elsif ($filename =~ /\,b?gz(ip)?$/) {
            my $zcat ||= first {-x $_} map {"$_/zcat"} @PATH;
            die("cannot find executable $zcat") if (! -x $zcat);
            my $fh_opener = "$zcat  $filename |";
            open $fh, $fh_opener;
        } else {
            open ($fh, "<$filename");
        }
    }
    return $fh;
}


sub check_GT_pedigree {
    my ( $progeny,$mother,$father ) = @_;
    my @mother_GT = sort (split /\/|\|/, $mother);
    my @father_GT = sort (split /\/|\|/,  $father);
    my @progeny_GT = sort (split /\/|\|/, $progeny);


    if ( ($mother =~ /\./ && $father =~ /\./) || $progeny =~ /\./ ) { #if both parents are not genotyped || progeny is not genotyped                                                                              
        return 0;
        last;
    } elsif ( $mother =~ /\./ || $father =~ /\./ ) {
        my $known_GT = first { $_ !~ /\./ } ($mother, $father);
        my %known_GT = map {$_ => 1} (split /\/|\|/, $known_GT);
        ( grep ($known_GT{$_}, @progeny_GT) ) ? return 1 : return 0;
    } else {
        my @combination_GT;
        foreach my $i (@mother_GT) {
            foreach my $j (@father_GT) {
                push @combination_GT, [sort ($i, $j)]
            }
        }
        ( grep  @{$_} ~~ @progeny_GT , @combination_GT ) ? return 1 : return 0;
#       return @progeny_GT;
    }
}



sub segCodeCheck {
    my ( $f1_1, $f1_2 ) = @_;
    my @f1_1_GT = sort (split /\/|\|/, $f1_1);
    my @f1_2_GT = sort (split /\/|\|/, $f1_2);
    my %f1_1_GT = map {$_ => 1} @f1_1_GT;
    my %f1_2_GT = map {$_ => 1} @f1_2_GT;

                                                                                                                                                                    
    my $usable =( (keys %f1_1_GT == 1) && (keys %f1_2_GT == 1) ) ? 0 : 1;
    # two F1 are homozygous, not informative 
    my %counts;
    if ($usable) {
        my $uniqueNo = grep !$counts{$_}++, (@f1_1_GT, @f1_2_GT);

        given ($uniqueNo) {
            when(4) {
                # case of abXcd
		my ($a, $b) = @f1_1_GT;
                my ($c, $d) = @f1_2_GT;
                return 'abXcd';
            }
	    when(3) {
                # case of efXeg
		my @intersect =  grep ( $f1_2_GT{$_}, @f1_1_GT );
                my $e = join "", @intersect;
		my $f = join "", ( grep {$_ ne  $e} @f1_1_GT );
                my $g = join "", ( grep {$_ ne  $e} @f1_2_GT );
                return 'efXeg';
            }
	    when(2) {
                # case of hkXhk, lmXll, nnXnp
		if ( (join "", @f1_1_GT) eq  (join "", @f1_2_GT) ) {
                    # case of hkXhk
                    my ($h,$k) = @f1_1_GT;
                    return 'hkXhk';
                } elsif ( keys %f1_2_GT == 1 ) {
                    # case of lmXll
                    my $l = $f1_2_GT[0];
                    my $m =  join "", ( grep {$_ ne  $l} @f1_1_GT );
                    return 'lmXll';
                } elsif ( keys %f1_1_GT ==1 ) {
                    # case of nnXnp
		    my $n = $f1_1_GT[0];
                    my $p =  join "", ( grep {$_ ne  $n} @f1_2_GT );
                    return 'nnXnp';
                } else {
                    return 'something else';
                } # end of if-elsif-else
	    }
            default { print 'Everything else' }
        }  # end of swithch case
    } else {
        return 'not informative';
    }# end of if 
}


sub segCodeConvert {
    my $autoCor = 1;
    my ( $f1_1, $f1_2, @f2 ) = @_;
    my @f1_1_GT = sort (split /\/|\|/, $f1_1);
    my @f1_2_GT = sort (split /\/|\|/, $f1_2);
    my %f1_1_GT = map {$_ => 1} @f1_1_GT;
    my %f1_2_GT = map {$_ => 1} @f1_2_GT;
   
    
    my $usable =( (keys %f1_1_GT == 1) && (keys %f1_2_GT == 1) ) ? 0 : 1;
    # two F1 are homozygous, not informative                                                                                                                                                                    
	my %counts;
	if ($usable) {
	    my @f2_candidates;
	    my @f2_types;
	    my @f2_GT = map { (my $s=$_) = s/\/|\|//; $_ } @f2;
	    my @fp_homo; # false positive homozygote due to allele dropout

#	    for my $i (@f2_GT) {
#		print "$i\n";
#	    }

	    my $uniqueNo = grep !$counts{$_}++, (@f1_1_GT, @f1_2_GT);
	      
	    given ($uniqueNo) {
		when(4) {
                # case of abXcd                                                                                                                                                                                  
		    my ($a, $b) = @f1_1_GT;
		    my ($c, $d) = @f1_2_GT;
		    print "<abXcd>\t";
		    @f2_candidates = ( [$a,$c], [$a,$d], [$b,$c], [$b,$d] );
		    @f2_types = ("ac","ad","bc","bd");

		    for my $i (@f2_GT) {
                        # write F2 genotype code 
			my @f2_GT_ind = split "", $i;
			if ( $i =~ /\./ ) {
			    print "--\t";
			} elsif ( unorderCom(\@f2_GT_ind,$f2_candidates[0]) ) {
                            #print "@f2_GT_ind:@{$f2_candidates[0]}:$f2_types[0]\t";
			    print "$f2_types[0]\t";
                        } elsif ( unorderCom(\@f2_GT_ind,$f2_candidates[1]) ) {
                            #print "@f2_GT_ind:@{$f2_candidates[1]}:$f2_types[1]\t";
			    print "$f2_types[1]\t";
                        } elsif ( unorderCom(\@f2_GT_ind,$f2_candidates[2])  ) { 
			   #print "@f2_GT_ind:@{$f2_candidates[2]}:$f2_types[2]\t";
			    print "$f2_types[2]\t";
			} elsif ( unorderCom(\@f2_GT_ind,$f2_candidates[3]) ) {
			    #print "@f2_GT_ind:@{$f2_candidates[3]}:$f2_types[3]\t";
			    print "$f2_types[3]\t";
                        }  else {
                            print "--\t";
			}
		    } # end of for loop 
		} # end of when 4

		when(3) {
                # case of efXeg
		    my @intersect =  grep ( $f1_2_GT{$_}, @f1_1_GT );
		    my $e = join "", @intersect;
		    my $f = join "", ( grep {$_ ne  $e} @f1_1_GT );
		    my $g = join "", ( grep {$_ ne  $e} @f1_2_GT );
		    print "<efXeg>\t";
		    @f2_candidates = ( [$e,$e], [$e,$f], [$e,$g],[$f,$g] );
		    @f2_types = ("ee", "ef", "eg","fg");
		    
		    for my $i (@f2_GT) {
                        # write F2 genotype code 
			my @f2_GT_ind = split "", $i;
			if ( $i =~ /\./ ) {
                            print "--\t";
                        } elsif ( unorderCom(\@f2_GT_ind,$f2_candidates[0]) ) {
                            #print "@f2_GT_ind:@{$f2_candidates[0]}:$f2_types[0]\t";
			    print "$f2_types[0]\t";
                        } elsif ( unorderCom(\@f2_GT_ind,$f2_candidates[1]) ) {
			    #print "@f2_GT_ind:@{$f2_candidates[1]}:$f2_types[1]\t";
			    print "$f2_types[1]\t";
                        } elsif ( unorderCom(\@f2_GT_ind,$f2_candidates[2])  ) {
			    #print "@f2_GT_ind:@{$f2_candidates[2]}:$f2_types[2]\t";
			    print "$f2_types[2]\t";
			} elsif ( unorderCom(\@f2_GT_ind,$f2_candidates[3])  ) {
			    #print "@f2_GT_ind:@{$f2_candidates[3]}:$f2_types[3]\t";
			    print "$f2_types[3]\t";
                        }  else {
                            #print "@f2_GT_ind:--\t";
			    print "--\t";
			}
		    } # end of for loop 		    
		} # end of when 3

		when(2) {
		    # case of hkXhk, lmXll, nnXnp                        
		    my ($h, $k, $l, $m, $n, $p);
		    if ( (join "", @f1_1_GT) eq  (join "", @f1_2_GT) ) {
			# case of hkXhk 
			($h, $k) = @f1_1_GT;
			print "<hkXhk>\t";
			@f2_candidates = ([$h,$h], [$h,$k], [$k,$k]);
			@f2_types = ("hh","hk","kk");
		    } elsif ( keys %f1_2_GT == 1 ) {
			# case of lmXll
			$l = $f1_2_GT[0];
			$m =  join "", ( grep {$_ ne  $l} @f1_1_GT );
			print "<lmXll>\t";
			@f2_candidates = ([$l,$l], [$l,$m]);
			@f2_types = ("ll", "lm");
			@fp_homo = ($m, $m);
		    } elsif ( keys %f1_1_GT ==1 ) {
			# case of nnXnp                   
			$n = $f1_1_GT[0];
			$p =  join "", ( grep {$_ ne  $n} @f1_2_GT );
			print "<nnXnp>\t";
			@f2_candidates = ([$n,$n], [$n,$p]);
			@f2_types = ("nn", "np");
			@fp_homo = ($p,$p);
#			print "nn:$n.$n\tnp:$n.$p\n";
			#print "@{$f2_candidates[0]}\n";
		    } else {
			print 'something else';
			next;
		    } # end of if-elsif-else
		    
		    for my $i (@f2_GT) {
			# write F2 genotype code
			my @f2_GT_ind = split "", $i;
			if ( $i =~ /\./ ) {
                            print "--\t";
			} elsif ( unorderCom(\@f2_GT_ind,$f2_candidates[0]) ) {
			    #print "@f2_GT_ind:@{$f2_candidates[0]}:$f2_types[0]\t";
			    print "$f2_types[0]\t";
			} elsif ( unorderCom(\@f2_GT_ind,$f2_candidates[1]) ) {
			    #print "@f2_GT_ind:@{$f2_candidates[1]}:$f2_types[1]\t";
			    print "$f2_types[1]\t";
			} elsif ( (defined $f2_candidates[2])  && unorderCom(\@f2_GT_ind,$f2_candidates[2])  ) { # in case of three candidate progeny genotypes
			    #print "@f2_GT_ind:@{$f2_candidates[2]}:$f2_types[2]\t";
			    print "$f2_types[2]\t";
			}  elsif ( ($autoCor == 1) &&  @fp_homo && unorderCom(\@f2_GT_ind,\@fp_homo) ) {
			    #print "@f2_GT_ind:@fp_homo:", uc $f2_types[1],"\t";
			    print uc $f2_types[1], "\t";
			}  else {
			    #print "@f2_GT_ind:--\t";
			    print "--\t";
			}
		    } # end of for loop

		} # end of when (2)

		default { print 'Everything else' }
	    }  # end of switch case 

	    print "\n"; # end of print row
	} else {
	    print "not informative\n";
	}# end of if     
}




sub unorderCom {
    my $arra_ref = shift;
    my $arrb_ref = shift;
 
    die "arrays have different size" if ( (scalar @{$arra_ref}) != (scalar @{$arrb_ref}) );
    # assume @arra and @arrb have the same size
    my %hasha = map {$_ => 1} @{$arra_ref};
    my %hashb = map {$_ => 1} @{$arrb_ref};
    my @b_contain_a = grep { $hashb{$_} } @{$arra_ref};
    my @a_contain_b = grep { $hasha{$_} } @{$arrb_ref};
    (scalar @b_contain_a == scalar @a_contain_b && scalar @a_contain_b == scalar @{$arra_ref} ) ? return 1 : return 0;
} 

my @arr1 = qw /0 0/;
my @arr2 = qw /1 1/;
my $result = unorderCom (\@arr1, \@arr2);
print "$result\n";



1;
