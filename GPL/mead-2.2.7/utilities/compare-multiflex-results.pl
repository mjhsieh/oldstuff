#! /usr/bin/perl -w

$usage = "Usage:
$0 [-thresh f] rootname1 rootname2
Compares the results (.pkint and .g files) of multiflex runs corresponding
to rootname1 and rootname2.  The rootnames are formed as, dirname/molname,
where dirname is the directory in which the run was done, and molname
is the molname argument given to multiflex.  The two runs must have
identical names, types, and orderings of the sites.  Only differences
greater than 0.01 pK units are reported.  The option, -thresh f, changes
the reporting threshold to f.
";

$ln10kT = 0.00413391;
$thresh = 0.01;

while (@ARGV) {
    $arg = shift @ARGV;
    if ($arg =~ /^-th/) {
	$thresh = shift @ARGV;
	die "Bad threshold value, $thresh\n $usage" unless $thresh > 0;
	next;
    }
    if (defined($rootname1) && defined($rootname2)) {
	warn "Excess command arguments ingorned.\n $usage";
	last;
    }
    elsif (defined($rootname1)) {
	$rootname2 = $arg;
    }
    else {
	$rootname1 = $arg;
    }
}

(defined($rootname1) && defined($rootname2))
    || die "Missing one or both rootnames.\n $usage";

$pkint1_fname = "$rootname1.pkint";
$pkint2_fname =  "$rootname2.pkint";
$g1_fname = "$rootname1.g";
$g2_fname =  "$rootname2.g";

# $pkint1_fname = "/home/dillet/WTSRI/thio/ftp/1xoa/xoa1/xoa1.pkint";
# $pkint2_fname =  "try/1xoa/xoa1/xoa1.pkint";
# $g1_fname = "/home/dillet/WTSRI/thio/ftp/1xoa/xoa1/xoa1.g";
# $g2_fname =  "try/1xoa/xoa1/xoa1.g";


open (PKINT1, $pkint1_fname) || die "Couldn't open $pkint1_fname\n $usage";
open (PKINT2, $pkint2_fname) || die "Couldn't open $pkint2_fname\n $usage";

while(defined($line1 = <PKINT1>)) {
    chop $line1;
    @line1 = split / +/, $line1;
    @line1 == 3 || die "Wrong number of fields in line, \"$line1\"\n \
                        from file, $pkint1_fname.\n";
    defined($line2 = <PKINT2>)
	|| die "Did not find in file, $pkint2_fname, a line \
                corresponding to \"$line1\" from $pkint1_fname\n";
    chop $line2;
    @line2 = split / +/, $line2;
    @line2 == 3 || die "Wrong number of fields in line, \"$line2\"\n \
                        from file, $pkint2_fname.\n";
    die "Mismatch in site info between lines, \"$line1\" from $pkint1_fname, \
         and \"$line2\" from $pkint1_fname.\n"
	     unless $line1[1] eq $line2[1] && $line1[2] eq $line2[2];
    push @sitename, $line1[2];
    $pkint_diff = abs($line1[0] - $line2[0]);
    if ($pkint_diff > $thresh) {
	print "$line1[2]    pKint diff =  $pkint_diff\n";
    }
}
close(PKINT1);
close(PKINT2);
    
open (G1, $g1_fname) || die "Couldn't open $g1_fname\n $usage";
open (G2, $g2_fname) || die "Couldn't open $g2_fname\n $usage";

while(defined($line1 = <G1>)) {
    chop $line1;
    @line1 = split / +/, $line1;
    @line1 == 3 || die "Wrong number of fields in line, \"$line1\"\n \
                        from file, $g1_fname.\n";
    defined($line2 = <G2>)
	|| die "Did not find in file, $g2_fname, a line \
                corresponding to \"$line1\" from $g1_fname\n";
    chop $line2;
    @line2 = split / +/, $line2;
    @line2 == 3 || die "Wrong number of fields in line, \"$line2\"\n \
                        from file, $g2_fname.\n";
    die "Mismatch in site info between lines, \"$line1\" from $g1_fname, \
         and \"$line2\" from $g1_fname.\n"
	     unless $line1[0] == $line2[0] && $line1[1] == $line2[1];
    $g_diff = abs($line1[2] - $line2[2])/$ln10kT;
    next if $line1[0] > $line1[1]; # Skip redundancies in symm matrix
    if ($g_diff > $thresh) {
	print "${sitename[$line1[0]]} ${sitename[$line1[1]]} "
	    . "interaction diff =   $g_diff pK units\n";
    }
}
close(G1);
close(G2);


