#!/usr/bin/perl -w

$chmapfile = "parse_charges";
$radmapfile = "parse_radii";

$usage = "Usage: $0 [-chmap cfile] [-radmap rfile] [infile.pdb] Reads
the .pdb named on the command line or pdb data from standard input,
and writes .pqr to standard output with charges and radii
assigned according to the files specified with the -chmap and -radmap
options.  If no chmap or radmap files are specified, the files,
$chmapfile and $radmapfile in the current directory are used.
These default file paths can easily be changed by editing the top
few lines of the script.
";

# Here put residue name translations to be applied to residue
# names in the input pdb file prior to lookup in the parse files.

%resnametrans = ("AMN", "BKN",
		 "CBX", "BKC",
                 "HSP", "HI+",
                 "HSC", "HI+");

%atnametrans = ("H", "HN",
		"CL", "CA");

%keytrans = ("CD:ILE", "CD1:ILE",
             "HG1:SER", "HG:SER",
             "HG1:CYS", "HG:CYS");

while (@ARGV > 1) {
    $w = shift @ARGV;
    if ($w eq "-radmap") {
	@ARGV || die "missing file name after -radmap option\n";
	$radmapfile = shift @ARGV;
    }
    elsif ($w eq "-chmap") {
	@ARGV || die "missing file name after -chmap option\n";
	$chmapfile = shift @ARGV;
    }
}



open (RADS, $radmapfile) || die "Could not open file $radmapfile\n$usage";
while ($_ = <RADS>) {
    s/!.*$//;
    s/\s+$//;
    split;
    next if @_ < 3;
    $radmap{"$_[0]:$_[1]"} = $_[2];
}

open (CHS, $chmapfile)  || die "Could not open file $chmapfile\n$usage";;
while ($_ = <CHS>) {
    s/!.*$//;
    s/\s+$//;
    split;
    next if @_ < 3;
    $chmap{"$_[0]:$_[1]"} = $_[2];
}

while (<>) { # PDB reading, PQR writing loop
    split;
    if (@_ >= 10 && $_[0] eq "ATOM") {
	$atname = $_[2];
	$resname = $_[3];
	$resname = $resnametrans{$resname} if defined $resnametrans{$resname};
	$atname = $atnametrans{$atname} if defined $atnametrans{$atname};
	$key = "$atname:$resname";
        $key = $keytrans{$key} if defined $keytrans{$key};
	if (defined($radmap{$key})) {
	    $radius = $radmap{$key};
	}
	else {
	    warn "No radius found for $key";
	    $radius = 0;
	}
	if (defined($chmap{$key})) {
	    $charge = $chmap{$key};
	}
	else {
	    warn "No charge found for $key";
	    $charge = 0;
	}

	printf "ATOM  %6d %4s %4s  %5d  %7.3f %7.3f %7.3f %6.3f %6.3f\n",
	    $_[1], $atname, $resname, $_[4], $_[5], $_[6], $_[7],
	    $charge, $radius;
    }
    else {
	print;
    }
}
