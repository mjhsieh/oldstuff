#!/usr/bin/perl

# Charge and radius sanity check for .pqr file.

$nocheckh = 1;

$first = 1;
while (<>) {
    /^ATOM/ || next;
    @Fld = split(' ');
    $resnum = $Fld[4];
    $resname = $Fld[3];
    $q = $Fld[8];
    $r = $Fld[9];
    if ($r > 5.0 || $r < 1.0) {
	if (!($nocheckh &&  $Fld[2] =~ /^H/)) {
	    print "WARNING: Weird radius:\n";
	    print;
	}
    }
    if ($resnum != $prevnum && !$first) {
	printf "%4d %4s charge = %7.3f    rmscharge = %7.3f\n",
	    $prevnum, $prevname, $qtot, sqrt($qsqtot/$n);
	$qtot = $qsqtot = $n = 0;
    }
    $qtot += $q;
    $qsqtot += $q * $q;
    ++$n;
    $prevnum = $resnum;
    $prevname = $resname;
    $first = 0;
}

printf "%4d %4s charge = %7.3f    rmscharge = %7.3f\n",
    $prevnum, $prevname, $qtot, sqrt($qsqtot/$n);
$qtot = $qsqtot = $n = 0;

