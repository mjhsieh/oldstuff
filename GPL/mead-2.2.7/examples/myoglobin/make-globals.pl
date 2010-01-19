#! /usr/bin/perl

# construct a global g-matrix and pkint file from its pieces
# Usage: make-globals <molecule-name>

# This script munges the pkint values and site-site interaction
# matrices from the multimead runs of the various tautomers into
# "global" pkint values and site-site interactions so that each
# histidine is represented by two psuedo-sites: one for each tautomer.
# It is necessary to put in some big interaction to suppress a "doubly
# deprotonated" state, and to adjust to other pkhalf values to stop
# sites from being seen twice.

# The resulting pkint values and site-site interactions (the .g file)
# should then be fed to some sort of multi-site titration solver, like
# Paul Berozas mcti.  The protonation values coming out of that solver
# for histidine will need to be de-munged to get protonation and
# tautomer populations at each pH.


#   nw = # of non-his sites:
#   nhis = # of his sites
#   input files are named $prefix-{del,eps,deln}.g etc.

$prefix=$ARGV[0];
if ($prefix =~ /MetMut/) {
  $nw = 50;
}
else {
  $nw = 51;
}
$nhis = 11;
$global_dim = $nw + $nhis + $nhis;

if (-e "$prefix.global.pkint" || -e "$prefix.global.g") {
	die "won't clobber existing files!\n";
}
$conv=241.919;
open (PKINT, ">$prefix.global.pkint");

#   read in (uncorrected) pk-intrinsic values:
open (IN,"${prefix}-del.pkint") || die "can't find ${prefix}-del.pkint\n";
while (<IN>) {
    $n++;
    ($pk[$n],$cat[$n],$name[$n]) = split(' ');
}
close IN;
open (IN,"${prefix}-eps.pkint") || die "can't find ${prefix}-eps.pkint\n";
while (<IN>) { 
    if ($. > $nw) {
	$n++;
	($pk[$n],$cat[$n],$name[$n]) = split(' ');
    }
    else {
	($pke[$.]) = split(' ');
    }
}
close IN;

#  read in the entire delta g-matrix making some pkint adjustments as we go.
open (IN,"${prefix}-del.g") || die "can't find ${prefix}-del.g\n";
while (<IN>) { 
    ($i,$j,$w) = split(' ');
    $g_global[($i-1)*$global_dim + $j-1] = $w;
    if ($j > $nw) {
	$pke[$i] += $w*$conv;
    }
}
close IN;

# pick out the relevant portions of the "epsilon" file wile adusting pkints
open (IN,"${prefix}-eps.g") || die "can't find ${prefix}-eps.g\n";
while (<IN>) {
	($i,$j,$w) = split(' ');
	if ($i <= $nw && $j >  $nw) { 
        	$g_global[($i-1)*$global_dim + $j+$nhis-1] = $w;
		$pk[$i] += $w*$conv;
	}
	if ($i > $nw && $j <= $nw) {
	    $g_global[($i+$nhis-1)*$global_dim + $j-1] = $w;
	}
	if ($i > $nw && $j > $nw) { 
        	$g_global[($i+$nhis-1)*$global_dim + $j+$nhis-1] = $w;
	}
}
close IN;

#  for each of the $prefix.g.del# files, read the appropriate files
#    and add to that block of the output file:

$wdiag = 20.0/$conv;		# A big value to forbid double deprotonation
for ($n=1; $n<=$nhis; $n++) {
    open (IN,"$prefix-del$n.g") || die "can't find $prefix-del$n.g\n";
    while (<IN>) {
	($i,$j,$w) = split(' ');
	if ($i==$n && $i != $j) {
	    $g_global[($i+$nw-1)*$global_dim + $j+$nw+$nhis-1] = $w;
	    $g_global[($j+$nw+$nhis-1)*$global_dim + $i+$nw-1] = $w;
	    $pk[$i+$nw] += $w*$conv;
	    $pk[$j+$nw+$nhis] += $w*$conv;
	}
    }
    close IN;

#  finally, add the large diagonal elements:
    $g_global[($n+$nw-1)*$global_dim + $n+$nw+$nhis-1] = $wdiag;
    $g_global[($n+$nw+$nhis-1)*$global_dim + $n+$nw-1] = $wdiag;
    $pk[$n+$nw] += $wdiag*$conv;
    $pk[$n+$nw+$nhis] += $wdiag*$conv;
}

# Print out the new g file
open (G,">$prefix.global.g");
for ($i = 1; $i <= $global_dim; ++$i) {
   for ($j = 1; $j <= $global_dim; ++$j) {
      printf G "%4d %4d %15.6e\n", $i, $j, $g_global[($i-1)*$global_dim + $j-1];
  }
}
close G;

#  print out the new pkint file:
for ($n=1;$n<=$nw+$nhis+$nhis;$n++) { 
	printf PKINT "%10.5f %1s %s\n",$pk[$n],$cat[$n],$name[$n];
}
close PKINT;
