package vcf2qtl;

use strict;
use warnings;
use PerlIO::gzip;
use Env qw ( @PATH );
use List::Util qw ( first max );
use Data::Dumper;
use feature qw(switch);
use Exporter qw(import);

our @EXPORT_OK = qw (open_input uniq phase_check phase transform first max switch @PATH);

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


sub uniq {
    my %seen;
    return grep { !$seen{$_}++ } @_;
}



sub phase_check {
    my ( $ref_progenies, $mother, $father) = @_;
    my @mother_GT = sort (split /\/|\|/, $mother);
    my @father_GT = sort (split /\/|\|/,  $father);
    foreach my $progeny ( @$ref_progenies ) {

	my @progeny_GT = sort (split /\/|\|/, $progeny);
	print "progeny:$progeny\tm:$mother\tf:$father\t";
	# Case 1: in case of being unable to phase
	# f1 homozygotes or two f0 have the same genotype
	if ( keys %{ { map{$_ ,1} @progeny_GT } } == 1 or join (".", @mother_GT) eq join (".", @father_GT) ) {
	    print "unphasable\t";
	} else { ## Case 2: Can be phased as f1 is heterozygous
	    my %pos_all; # assume it is diploid
	    $pos_all{'0'} = $progeny_GT[0];
	    $pos_all{'1'} = $progeny_GT[1];

	    # go through each allele
	    foreach my $key (sort {$a <=> $b} keys %pos_all) { # position 0 and position 1
		# In case 2, allele 0 from mother
		if ( (grep /$pos_all{$key}/, @mother_GT) && !(grep /$pos_all{$key}/, @father_GT) ) {
		    $pos_all{'0'} = $pos_all{$key};
		    $pos_all{'1'} = join "", grep {$_ ne $pos_all{$key}} @progeny_GT;
		    ( @progeny_GT ~~ @mother_GT && $father eq "." ) ? print "unphasable:$progeny\t" : 
			print "phased_progeny:", join ("/", ($pos_all{'0'},$pos_all{'1'})), "\t";
		    last;
		    # In case 2, allele 0 from father
		} elsif ( (grep /$pos_all{$key}/, @father_GT) && ! (grep /$pos_all{$key}/, @mother_GT) ) {
		    $pos_all{'1'} = $pos_all{$key};
		    $pos_all{'0'} = join "", grep {$_ ne $pos_all{$key}} @progeny_GT;
		    ( @progeny_GT ~~ @father_GT && $mother eq "." ) ? print "unphasable:$progeny\t" : 
			print "phased_progeny:", join ("/", ($pos_all{'0'},$pos_all{'1'})), "\t";
		    last;
		} else {
		    print "violate\t" if ($key eq '1');
		    next;
		} # end of if
	    } # end of foreach
	    
	} # end of if
	print "\n"; # go to next individual
    } # end of progenies loop
}

sub phase {
    my ( $ref_progenies, $mother, $father) = @_;
    my @mother_GT = sort (split /\/|\|/, $mother);
    my @father_GT = sort (split /\/|\|/,  $father);
    my @phased_progenies;
    my $phased_progeny;
    
    foreach my $progeny ( @$ref_progenies ) {
        
        my @progeny_GT = sort (split /\/|\|/, $progeny);
        
        # Case 1: in case of being unable to phase as
        # progeny's genotype is homozygotes or missing, but is actually phased already
        if ( keys %{ { map{$_ ,1} @progeny_GT } } == 1 ) {
            
            #print "unphasable\t";
            $phased_progeny = join( "|",@progeny_GT);
            
            # case 1: in case of being unable to phase as
            # two parents' have the same genotype
        } elsif ( join (".", @mother_GT) eq join (".", @father_GT) ) {
            
            #print "unphasable\t";
            $phased_progeny = $progeny;
            
        } else { ## Case 2: Can be phased as progeny is heterozygous
            
            my %pos_all; # assume it is diploid
            $pos_all{'0'} = $progeny_GT[0];
            $pos_all{'1'} = $progeny_GT[1];
            
            # go through each allele
            foreach my $key (sort {$a <=> $b} keys %pos_all) { # position 0 and position 1
                # In case 2, allele 0 from mother
                if ( (grep /$pos_all{$key}/, @mother_GT) && !(grep /$pos_all{$key}/, @father_GT) ) {
                    $pos_all{'0'} = $pos_all{$key};
                    $pos_all{'1'} = join "", grep {$_ ne $pos_all{$key}} @progeny_GT;
                    # Case 3: progeny is the same of one of f1, the other f1 genotype is unknown
                    ( (@progeny_GT ~~ @mother_GT) && ($father =~ /^\./) ) ?
                    ( $phased_progeny =  $progeny ):
                    ( $phased_progeny = join ( "|", ($pos_all{'0'},$pos_all{'1'}) ) );
                    last;
                    # In case 2, allele 0 from father
                } elsif ( (grep /$pos_all{$key}/, @father_GT) && ! (grep /$pos_all{$key}/, @mother_GT) ) {
                    $pos_all{'1'} = $pos_all{$key};
                    $pos_all{'0'} = join "", grep {$_ ne $pos_all{$key}} @progeny_GT;
                    # Case 3: progeny is the same of one of f1, the other f1 genotype is unknown
                    ( (@progeny_GT ~~ @father_GT) && ($mother =~ /^\./)  ) ?
                    ( $phased_progeny =  $progeny ):
                    ( $phased_progeny = join ( "|", ($pos_all{'0'},$pos_all{'1'}) )  );
                    
                    last;
                } else {
                    $phased_progeny = $progeny if ($key eq '1');
                    next;
                } # end of if-elsif-else
            } # end of foreach
            
        } # end of if-else
        push @phased_progenies, $phased_progeny;
        
    } # end of progenies loop
    return @phased_progenies;
    
}



sub transform {
    my ($ref_progenies, $mother, $father) = @_;
    my @mother_GT = split (/\/|\|/, $mother);
    my @father_GT = split (/\/|\|/, $father);
    my %fourwaycross; # it doesn't matter if one or two parents cannot be phased.
    my @transform_pro_geno;
    $fourwaycross{'A'} = $mother_GT[0];
    $fourwaycross{'B'} = $mother_GT[1];
    $fourwaycross{'C'} = $father_GT[0];
    $fourwaycross{'D'} = $father_GT[1];
    
    foreach my $progeny (@$ref_progenies) {
        (my $progeny_GT = $progeny) =~  s/\/|\|//g;
        my @progeny_GT = split ("", $progeny_GT);
        
        if ( $progeny =~ /\|/ || keys %{ { map{$_ ,1} @progeny_GT } } == 1 ) {# if progeny is phasable or (progeny is homozygous, actually considered in the phase function)
            
            if ( $progeny_GT =~ /^\./ ) {
                push @transform_pro_geno, "NA";
            } elsif (  $progeny_GT eq $fourwaycross{'A'} . $fourwaycross{'D'} && $progeny_GT eq  $fourwaycross{'B'} . $fourwaycross{'C'} ) {
                push @transform_pro_geno, "10";
            } elsif  ($progeny_GT eq $fourwaycross{'A'} . $fourwaycross{'C'} && $progeny_GT eq $fourwaycross{'B'} . $fourwaycross{'D'}) {
                push @transform_pro_geno, "9";
            } elsif ($progeny_GT eq $fourwaycross{'A'} . $fourwaycross{'D'} && $progeny_GT eq $fourwaycross{'B'} . $fourwaycross{'D'}) {
                push @transform_pro_geno, "8";
            } elsif ($progeny_GT eq $fourwaycross{'A'} . $fourwaycross{'C'} && $progeny_GT eq $fourwaycross{'B'} . $fourwaycross{'C'}) {
                push @transform_pro_geno, "7";
            } elsif ($progeny_GT eq $fourwaycross{'B'} . $fourwaycross{'C'} && $progeny_GT eq $fourwaycross{'B'} . $fourwaycross{'D'}) {
                push @transform_pro_geno, "6";
            } elsif ($progeny_GT eq $fourwaycross{'A'} . $fourwaycross{'C'} && $progeny_GT eq $fourwaycross{'A'} . $fourwaycross{'D'}) {
                push @transform_pro_geno, "5";
            } elsif ($progeny_GT eq $fourwaycross{'B'} . $fourwaycross{'D'}) {
                push @transform_pro_geno, "4";
            } elsif ($progeny_GT eq $fourwaycross{'A'} . $fourwaycross{'D'}) {
                push @transform_pro_geno, "3";
            } elsif ($progeny_GT eq $fourwaycross{'B'} . $fourwaycross{'C'}) {
                push @transform_pro_geno, "2";
            } elsif ($progeny_GT eq $fourwaycross{'A'} . $fourwaycross{'C'}) {
                push @transform_pro_geno, "1";
                #	} elsif ($progeny_GT ne $fourwaycross{'A'} . $fourwaycross{'C'}) {
                #	    push @transform_pro_geno, "11";
            } else {
                push @transform_pro_geno,"unknown";
            }   # end of if-elsif
            #else {
            #   my @progeny_GT = split ("", $progeny_GT);
            #   if ( ($progeny_GT[0] eq $progeny_GT[1]) &&
            
            #       push @transform_pro_geno, "unknown\t"
            
        } elsif ($progeny =~ /\// ) { # in case of progeny is not phsable
            
            if ( $progeny_GT =~ /^\./ ) {
                push @transform_pro_geno, "NA";
            } elsif ( ($fourwaycross{'A'} . $fourwaycross{'D'}) eq ($fourwaycross{'C'} . $fourwaycross{'B'}) && ($fourwaycross{'A'}, $fourwaycross{'D'}) ~~ @progeny_GT ) {
                push @transform_pro_geno, "10";
            } elsif ( ($fourwaycross{'A'} . $fourwaycross{'C'}) eq ($fourwaycross{'D'} . $fourwaycross{'B'}) && ($fourwaycross{'A'}, $fourwaycross{'C'}) ~~ @progeny_GT )  {
                push @transform_pro_geno, "9";
                #   } elsif ($progeny_GT ne $fourwaycross{'A'} . $fourwaycross{'C'}) {
                #       push @transform_pro_geno, "11";
            } else {
                push @transform_pro_geno, "na";
            }# end of if-elsif-else
            
        }# end of outer if
        
    } # end of foreach
    return @transform_pro_geno;
} # end of program


1;
