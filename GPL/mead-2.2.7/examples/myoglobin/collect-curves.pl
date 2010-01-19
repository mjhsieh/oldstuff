#! /usr/bin/perl  -w
#
#  read output from Monte Carlo titration program, print out
#     titration behavior in nice table on standard output
#
#  Usage:  pkhalf.pl <mcti-output-file>
#

# This is a descendent of the pkhalf and hisplot scripts by Dave Case.
# It his been simplified some and it is modified for the new output
# format of Berozas mcti program.  It puts all output in one file (stdout)
# rather than making separate hisplot files, and it no longer uses plot77 format.
# Eliminated width calculations.

$infile = $ARGV[0];

open (IN,"$infile") || die "Input file $infile not found\n";

# first line of input should have number of sites and pH points
$_ = <IN>;
($global_dim, $npH) = split (' ');
# Sanity checks
$global_dim > 0 || die "ERROR number of sites = $global_dim\n";
$npH > 0 || die "npH = $npH\n";

# Next come a listing of sites.

for ($i=0; $i<$global_dim; ++$i) {
    $_ = <IN> || die "Unexpected end of input\n";
    ($pkint, $cat, $name[$i]) = split (' ');
    $nameindex{$name[$i]} = $i;
}

# For each pH, read the pH value ... 
for ($ipH=0; $ipH<$npH; ++$ipH) {
    $_ = <IN> || die "Unexpected end of input\n";
    ($pH[$ipH]) = split (' ');
    $_ = <IN> || die "Unexpected end of input\n";
    ($jsite, $totp, $err) = split (' ');

# ... and then the protonation for each site.
    for ($isite=0; $isite<$global_dim; ++$isite) {
	$_ = <IN> || die "Unexpected end of input\n";
	($jsite, $prot[$ipH*$global_dim + $isite], $err) = split (' ');
	$jsite == $isite + 1 || warn "WARNING: Site number mismatch\n $isite\n $_";
    }
}

close IN;

#  now cycle through each residue and print the protonation state
#   at each pH found:

# Do the non-histidines first.

for ($isite=0; $isite<$global_dim; ++$isite) {
    if ($name[$isite] =~ /HIS/) {next;}
    print "------------------------------------------------------\n";
    print "Titration site: ${name[$isite]} \n";
    print "pH       prot.\n";
    $lastp = -1;
    $nhalf = 0;
    undef $pkhalf;
# loop through pH values
    for ($ipH=0; $ipH<$npH; ++$ipH) {
	$thisprot = $prot[$ipH*$global_dim + $isite];
	$thispH = $pH[$ipH];

	if ($thisprot < 0.98 && $thisprot > 0.02) {
	    printf "%7.2f   %7.4f\n", $thispH, $thisprot;
	}

#    two-point linear interpolation to get pkhalf:
	if ($lastp > 0.5 && $thisprot <= 0.5) {
	    $a = $lastp - 0.5; $b = 0.5 - $thisprot;
	    $pkhalf = ($a/($a+$b))*($thispH - $lastph) + $lastph;
	    $nhalf = $ipH;
	}
	$lastp = $thisprot;
	$lastph = $thispH;
    }

# Use six points to get Hill slope:
    undef $slope;
    if ($nhalf >= 4 ) {
	$sx = $sy = $sxx = $sxy = 0.0;
	for ($i = $nhalf-3; $i <= $nhalf +2; $i++) {
	    last if $i >= $npH;
	    $thisprot = $prot[$i*$global_dim + $isite];
	    $thispH = $pH[$i];
	    $y = log($thisprot / (1.0 - $thisprot))/2.303;
	    $sx += $thispH; $sxx += $thispH**2;
	    $sy += $y; $sxy += $y*$thispH;
	}
	$slope = -(6.0 * $sxy - $sx*$sy)/(6.0 * $sxx - $sx*$sx);
    }
    if (defined $pkhalf) {
	printf "pK(1/2) for %-11s  =  %7.3f\n", $name[$isite], $pkhalf;
    }
    if (defined $slope) {
	printf "Slope of Hill plot for %-11s  =  %7.3f\n", $name[$isite], $slope;
    }
}

# Loop through His and make hisplot files
foreach $delname (grep(/HISdel/, @name)) {
    $epsname = $delname;
    $epsname =~ s/del/eps/;
    $plainname = $delname;
    $plainname =~ s/del//;
    print "------------------------------------------------------\n";
    print "Titration site: $plainname\n";
    print "pH        prot.    del-taut.  eps-taut.\n";
    $idel = $nameindex{$delname};
    $ieps = $nameindex{$epsname};
    undef $pkhalf;
    $nhalf = 0;
    $lastp = -1;
    for ($ipH=0; $ipH<$npH; ++$ipH) {
# Here is where we un-munge the protonations of those funny double sites
# to get the tautomers.
	$del = 1.0 - $prot[$ipH*$global_dim + $idel];
	$eps = 1.0 - $prot[$ipH*$global_dim + $ieps];
	$thisprot = 1.0 - ($del + $eps);
	$thispH = $pH[$ipH];
	if ($thisprot < 0.98 && $thisprot > 0.02) {
	    printf "%7.2f   %7.4f   %7.4f   %7.4f\n", $thispH, $thisprot, $del, $eps;
	}

# Two-point liner interpolation as before.
	if ($lastp > 0.5 && $thisprot <= 0.5) {
	    $a = $lastp - 0.5; $b = 0.5 - $thisprot;
	    $pkhalf = ($a/($a+$b))*($thispH - $lastph) + $lastph;
	    $nhalf = $ipH;
	}
	$lastp = $thisprot;
	$lastph = $thispH;
    }
# end of pH loop

# Use six points to get Hill slope:
    undef $slope;
    if ($nhalf >= 4 ) {
	$sx = $sy = $sxx = $sxy = 0.0;
	for ($i = $nhalf-3; $i <= $nhalf +2; $i++) {
	    last if $i >= $npH;
	    $del = 1.0 - $prot[$i*$global_dim + $idel];
	    $eps = 1.0 - $prot[$i*$global_dim + $ieps];
	    $thisprot = 1.0 - ($del + $eps);
	    $thispH = $pH[$i];
	    $y = log($thisprot / (1.0 - $thisprot))/2.303;
	    $sx += $thispH; $sxx += $thispH**2;
	    $sy += $y; $sxy += $y*$thispH;
	}
	$slope = -(6.0 * $sxy - $sx*$sy)/(6.0 * $sxx - $sx*$sx);
    }
    if (defined $pkhalf) {
	printf "pK(1/2) for %-11s  =  %7.3f\n", $plainname, $pkhalf;
    }
    if (defined $slope) {
	printf "Slope of Hill plot for %-11s  =  %7.3f\n", $plainname, $slope;
    }
}
