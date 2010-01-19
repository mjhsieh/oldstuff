#!/usr/bin/perl -w

$usage = "Usage:
$0 -pdb pbd_file -psf psf_file -par par_file -top top_file [-hminrad f]

Reads pdb_file for coordiantes, a CHARMM psf_file for charges and atom
types, CHARMM param_file for radii, and CHARMM top_file for atom type
indices; then writes a .pqr file with charges and radii assigned to
standard output.  Optional flag, -hminrad f, makes f the minimum hydrogen
 radius.
";



@ARGV >= 8 || die $usage;

while(@ARGV) {
    $opt = shift @ARGV;
    if ($opt eq "-pdb") {
	$pdb_file = shift @ARGV;
	open (PDB, $pdb_file) || die "PDB file, $pdb_file not found\n\n$usage";
    }
    elsif ($opt eq "-psf") {
	$psf_file = shift @ARGV;
	open (PSF, $psf_file) || die "PSF file, $psf_file not found\n\n$usage";
    }
    elsif ($opt eq "-par") {
	$par_file = shift @ARGV;
	open (PAR, $par_file) || die "PAR file, $par_file not found\n\n$usage";
    }
    elsif ($opt eq "-top") {
	$top_file = shift @ARGV;
	open (TOP, $top_file) || die "TOP file, $top_file not found\n\n$usage";
    }
    elsif ($opt eq "-hminrad") {
	$hminrad = shift @ARGV;
	open (TOP, $top_file) || die "TOP file, $top_file not found\n\n$usage";
    }
    else {
	die "ERROR on command line near \"$opt\"\n\n$usage";
    }
}

die "ERROR: pdb file not specified\n\n$usage" unless defined ($pdb_file);
die "ERROR: psf file not specified\n\n$usage" unless defined ($psf_file);
die "ERROR: par file not specified\n\n$usage" unless defined ($par_file);
die "ERROR: top file not specified\n\n$usage" unless defined ($top_file);



warn "Scanning top file\n";

# We want the "MASS" entries which give the index for each atom type.
# Assume each index is on a single line of the form:
# MASS  index atomtype mass  [extra_junk]

while ($_ = <TOP>) {
    last if /^ *RESI/i;  # Assume that once RESI's start, no more MASSs
    next unless /^ *MASS/i;
    @_ = split;
    shift if $_[0] eq "";
    $atype_idx{$_[2]} = $_[1];
}
close (TOP);

# We want the radius for each atom type name.  These are in the
# NONBOND section lies between the NBON or NONB line and the NBFIX or HBON line

$nbon = 0;
while ($_ = &charmmline (PAR)) {
    /^ *(NONB|NBON)/i && do {$nbon=1; next;};
    /^ *NBFIX/i && do {$nbon = 0;};
    /^ *HBOND/i && do {$nbon = 0;};
    /^ *END/i && do {$nbon = 0;};
    next unless $nbon;
    @_ = split;
    shift @_ if $_[0] eq "";
    next if @_ < 4;
    $saw_matches=0;
    if ($_[0] =~ /[\*%\#\+]/) {
	# The atom spec is a CHARMM regex, so must translate it to perl regex.
	$preg = "^";
	for ($i=0; $i<length($_[0]); ++$i) {
	    $c = substr($_[0], $i, 1);
	  PRCHR: {
	      $preg .= "\\w*", last PRCHR if $c eq "*";
	      $preg .= "\\w", last PRCHR if $c eq "%";
	      $preg .= "\\d*", last PRCHR if $c eq "#";
	      $preg .= "\\d", last PRCHR if $c eq "+";
	      $preg .= $c;
	  }
	}	   
	$preg .= "\$";
	foreach $atype (keys %atype_idx) {
	    next unless $atype =~ /$preg/i;
	    if (defined($radlist[$atype_idx{$atype}])) {
		warn "Radius for $atype, which matches $_[0], already defined\n";
	    }
		
	    $radlist[$atype_idx{$atype}] = $_[3];
            $saw_matches = 1;
        }	    
	$saw_matches
	    || warn "WARNING: Nothing matches param NBON entry, $_[0]\n";
    }
    else {
	if (defined ($atype_idx{$_[0]})) {
	    $radlist[$atype_idx{$_[0]}] = $_[3];
	    $saw_matches = 1;
	}
    }
    if ($saw_matches==0) {
	warn "WARNING:\n";
	warn "Atom type, $_[0], seen in a param file NONBOND entry,\n";
	warn "does not match any type seen in top file MASS entries.\n\n";
    }
}

close (PAR);

$a=$b=0;  # supress stupid warning
sub byatidx { $atype_idx{$a} <=> $atype_idx{$b}; }
foreach $atype (sort byatidx keys %atype_idx) {
    $idx = $atype_idx{$atype};
    defined($radlist[$idx])
	|| do {warn "No radius found for atom type, $atype\n"; next;};
    if (defined ($hminrad) && $atype =~ /^H/i && $radlist[$idx] < $hminrad) {
	$radlist[$idx] = $hminrad;
    }
    warn "$atype  $atype_idx{$atype}  $radlist[$atype_idx{$atype}]\n";
}

# The atom entries of the PSF starts after the second occurence of
# an integer on a line by itself (after comment stripping, etc.).
# This integer will be the number of atom lines to follow.

$intseen = 0;

while ($_ = &charmmline(PSF)) {
    ++$intseen if /^\s*(\d+)\s*$/;
    ($natom = $1, last) if $intseen == 2;
}

warn "natom = $natom\n";

for ($i=0; $i<$natom; ++$i) {
    ($_ = <PSF>) || die "Unexpected EOF in PSF\n";
    @_ = split;
    shift @_ if $_[0] eq "";
    $_[0] == $i+1 || warn "Atom number mismatch, $_[0] != $i + 1\n";
    $charge[$i] = $_[6];
    $idx = $_[5];
    if (defined ($radlist[$idx])) {
	$radius[$i] = $radlist[$idx];
    }
    else {
	warn "WARNING: Unknown atom index, $idx, seen in PSF line,\n$_";
	warn "         Assigning zero radius to it.\n";
	$radius[$i] = 0;
    }
}

while ($_ = <PDB>) {
    @_ = split;
    if (@_ >= 10 && $_[0] eq "ATOM") {
	$i = $_[1] - 1;
	$atname = $_[2];
	$resname = $_[3];
	printf "ATOM  %6d %4s %4s  %5d  %7.3f %7.3f %7.3f %6.3f %6.3f\n",
	    $_[1], $atname, $resname, $_[4], $_[5], $_[6], $_[7],
	    $charge[$i], $radius[$i];
    }
    else {
	print;
    }
}

exit;

sub charmmline {
    local ($handle, $line);
    $handle = pop @_;
    do {
	($line = <$handle>) || return "";
	$line =~ s/!.*$//;
	while ($line =~ /- *$/) {
	    $line =~ s/- *$/ /;
	    chop $line;
	    ($line .= <$handle>) ||
		warn "WARNING: charmmline gets unexpected EOF from $handle\n";
	    $line =~ s/!.*$//;
	}
    } until $line =~ /\S/;
    $line;
}
    

