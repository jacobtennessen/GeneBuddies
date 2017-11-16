#!/usr/bin/perl -w
use strict;
use warnings;
use Getopt::Std;

##### ##### ##### ##### #####

use vars qw($opt_o $opt_s $opt_n $opt_r $opt_g $opt_c $opt_m);

# Usage
my $usage = "
BalanceLinkage.pl - simulates two linked polymorphisms under balancing selection with no epistasis
Copyright (C) 2017 by Jacob A Tennessen

This program is free software: you can redistribute it and/or modify
it under the terms of the GNU General Public License as published by
the Free Software Foundation, either version 3 of the License, or
(at your option) any later version.

This program is distributed in the hope that it will be useful,
but WITHOUT ANY WARRANTY; without even the implied warranty of
MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
GNU General Public License for more details.

You should have received a copy of the GNU General Public License
along with this program.  If not, see <http://www.gnu.org/licenses/>.

Usage: perl BalanceLinkage.pl options
 optional:
  -o    output file for R script to make plot of allele frequency change
  -s	selection coefficient (advantage to heterozygotes) [default = 0.005]
  -n    number of diploid individuals (Ne = 2N is twice this number) [default = 1000]
  -r    population-scaled recombination rate (rho = 4Ner) [default = 0]
  -g    number of generations to simulate [default = 100000]
  -c    number of replicate simulations [default = 1000]
";

#############

# command line processing.
getopts('o:s:n:r:g:c:');

my ($outputfile, $selection, $basepopsize, $rho, $gens, $reps);

if (defined $opt_o) {
    $outputfile = $opt_o;
}
if (defined $opt_s) {
    $selection = $opt_s;
} else {
    $selection = 0.005; 
}
if (defined $opt_n) {
    $basepopsize = $opt_n;
} else {
    $basepopsize = 1000;
}
if (defined $opt_r) {
    $rho = $opt_r;
} else {
    $rho = 0;
}
if (defined $opt_g) {
    $gens = $opt_g;
} else {
    $gens = 100000;
}
if (defined $opt_c) {
    $reps = $opt_c;
} else {
    $reps = 1000;
}

my $recrate = ($rho/(2*$basepopsize))/4;

print "Examining $reps replicates of $gens generations and $basepopsize diploid individuals.\n";
print "Selection coefficient is $selection.\n";  
print "Recombination rate is $recrate\n";

my $fixed = 0;
my $opposites = 0;
my $similarpair = 0;
my $all = 0;
my $three = 0;

my @rout;

my @firstfixed;

push @rout, "plot(-1000,-1000,xlim=c(0,$gens),ylim=c(0,1),xlab=\"Generations\",ylab=\"Frequency\")";

for (my $r = 1; $r <= $reps; $r++) {
    
    #print "\tRep $r\n";
    
    my @chartfreqs;
        
    my @gens;
    
    my %freqs;
    
    my $firstfixed;
    
    $freqs{"0\t0"} = 0.25;
    $freqs{"1\t1"} = 0.25;
    $freqs{"0\t1"} = 0.25;
    $freqs{"1\t0"} = 0.25;
    
    for (my $g = 0; $g <= $gens; $g++) {
        
        my $roundedg = (int($g/100))*100;
        
        my $freqfirst = $freqs{"0\t0"} + $freqs{"0\t1"};
        
        my $freqsecond = $freqs{"0\t0"} + $freqs{"1\t0"};
        
        if (($freqfirst == 1)||($freqfirst == 0)) {
            unless (defined $firstfixed) {
                $firstfixed = $g;
            }
            if (($freqsecond == 1)||($freqsecond == 0)) {
                #print "Ending at $g\n";
                last;
            }
        }
        
        if ((defined $outputfile)&&($roundedg == $g)) {
            push @chartfreqs, sprintf "%.4f", $freqfirst;         
            push @gens, $g;
        }
        
        my %gametes = (
            "0\t0" => 0,
            "0\t1" => 0,
            "1\t0" => 0,
            "1\t1" => 0,
        );
        
        for (my $i = 0; $i <= $basepopsize; $i++) {
            my $alleleA;
            my $alleleB;
            my $alleleArand = rand 1;
            if ($alleleArand <= $freqs{"0\t0"}) {
                $alleleA = "0\t0";
            } elsif ($alleleArand <= ($freqs{"0\t0"}+$freqs{"0\t1"})) {
                $alleleA = "0\t1";
            } elsif ($alleleArand <= ($freqs{"0\t0"}+$freqs{"0\t1"}+$freqs{"1\t0"})) {
                $alleleA = "1\t0";
            } else {
                 $alleleA = "1\t1";
            }
            my $alleleBrand = rand 1;
            if ($alleleBrand <= $freqs{"0\t0"}) {
                $alleleB = "0\t0";
            } elsif ($alleleBrand <= ($freqs{"0\t0"}+$freqs{"0\t1"})) {
                $alleleB = "0\t1";
            } elsif ($alleleBrand <= ($freqs{"0\t0"}+$freqs{"0\t1"}+$freqs{"1\t0"})) {
                $alleleB = "1\t0";
            } else {
                 $alleleB = "1\t1";
            }
            if ($alleleA =~ /^$alleleB$/) {
                $gametes{$alleleA} += 2;
            } else {
                my @allelea = split "\t", $alleleA;
                my @alleleb = split "\t", $alleleB;
                if (($allelea[0] == $alleleb[0])||($allelea[1] == $alleleb[1])) {
                    $gametes{$alleleA} += 1 + $selection;
                    $gametes{$alleleB} += 1 + $selection;
                } else {
                    $gametes{$alleleA} += (1 + ($selection*2))*(1-$recrate);
                    $gametes{$alleleB} += (1 + ($selection*2))*(1-$recrate);
                    my $altA = "$allelea[0]\t$alleleb[1]";
                    my $altB = "$alleleb[0]\t$allelea[1]";
                    $gametes{$altA} += (1 + ($selection*2))*$recrate;
                    $gametes{$altB} += (1 + ($selection*2))*$recrate;
                }
            }
        }
        
        my $popsize = $gametes{"0\t0"} + $gametes{"0\t1"} + $gametes{"1\t0"} + $gametes{"1\t1"};
        
        $freqs{"0\t0"} = $gametes{"0\t0"}/$popsize;
        $freqs{"0\t1"} = $gametes{"0\t1"}/$popsize;
        $freqs{"1\t0"} = $gametes{"1\t0"}/$popsize;
        $freqs{"1\t1"} = $gametes{"1\t1"}/$popsize;
        
    }
    
    if (defined $firstfixed) {
        push @firstfixed, $firstfixed;
    }
    
    if (defined $outputfile) {
        my $chartfreqs = join ",", @chartfreqs;
        my $gens = join ",", @gens;  
        push @rout, "lines(c($gens),c($chartfreqs),type=\"l\")";
    }
    
    my @retained;
    
    foreach my $f (keys %freqs) {
        unless ($freqs{$f} == 0) {
            push @retained, $f;
        }
    }
    
    if ((scalar(@retained)) == 1) {
        $fixed +=1;
    } elsif ((scalar(@retained)) == 3) {
        $three +=1;
    } elsif ((scalar(@retained)) == 4) {
        $all +=1;
    } else {
        my @firstallele = split "\t", $retained[0];
        my @secondallele = split "\t", $retained[1];
        if (($firstallele[0] == $secondallele[0])||($firstallele[1] == $secondallele[1])) {
            $similarpair +=1;
        } else {
            $opposites +=1;
        }
    }
}

print "Fixed: $fixed\nSimilar Two: $similarpair\nOpposite Two: $opposites\nThree: $three\nAll: $all\n";

my $middlepoint = $reps/2;

if ((scalar(@firstfixed)) >= $middlepoint) {
    @firstfixed = sort by_number @firstfixed;
    my $roundedmiddle = int($middlepoint-0.1);
    my $median = $firstfixed[$roundedmiddle];
    print "Median fixation time for first locus is $median generations.\n";
} else {
    print "Cannot calculate median fixation time: over half of simulations did not fix.\n";
}

if (defined $outputfile) {
    my $rout = join "\n", @rout;
    
    unless ( open(OUT, ">$outputfile") ) {
        print "Cannot open file \"$outputfile\" to write to!!\n\n";
        exit;
    }
    print OUT $rout;
    
    close (OUT);
}

####################

sub by_number {
    if ($a < $b) {-1} elsif ($a > $b) {1} else {0}
}
