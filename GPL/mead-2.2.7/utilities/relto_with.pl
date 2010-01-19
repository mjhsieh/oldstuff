#!/usr/bin/perl -w

# (c) Don Bashford, 

# This program allows you to find the energy (in pK units) that it
# takes to change the protonation microstate of selected groups of a
# multisite molecule if other sites were held in a specified
# protonation state.

# Usage: relto_with.pl molname -of [-p querysite] [-u querysite] ..
#                           -relto [-p querysite] [-u querysite] ..
#                           -with [-p fixsite ...] [-u fixsite ...]

# where querysite are names of a sites to be reported on and fixsite
# is the name of the site whose protonation state is to be fixed
# as protonated or unprotonated as specified by the -p or -u flags,
# respectively.  

# It has been creating by recuction from a more complex version that
# was specialized to the myolgobin calculations in which HIS had two
# neutral tautomers.


# A regular expression can be used for 
# fixsite so that groups of residues (like all ASPs can be easily
# specified.

# For sites not specified as fixsites, the protonation state
# corresponding to the "background" assumed by multimead (i.e.,
# sites in their neutral state) is used. 

# The program gathers data from mutlimead .pkint and .g files.

$prefix = $ARGV[0];  # A.K.A. molname
$pkint_file = "$prefix.pkint";
$g_file = "$prefix.g";
$conv = 241.902;


open (PKINT, $pkint_file) || die "Couldn't open $pkint_file\n";
$numsites = 0;
while(<PKINT>) {
    chop;
    ($pkint[$numsites], $c_or_a[$numsites], $nm) = split(" ");
    $nameindex{$nm} = $numsites;
    $sitename[$numsites] = $nm;
    $groupname[$numsites] = "";
    ++$numsites;
}
close PKINT;



# Read in the site-site interaction matrix, or G matrix

open (GF, $g_file) || die "Couldn't open $g_file\n";
while (<GF>) {
  ($i, $j, $f) = split (" ");
  $g[($i-1)*$numsites + ($j-1)] = $f;   #  2D arrays, C-style (yuk)
}


#parse the rest of the command line

# process the list of observed sites after the "-of" flag

$ARGV[1] eq "-of" || die "-of flag expected after $prefix, but not seen\n";

for ($argn = 2; $argn<@ARGV && $ARGV[$argn] ne "-relto" ; $argn += 2) {
  if ($ARGV[$argn] eq "-p") {$st = "p";}
  elsif ($ARGV[$argn] eq "-u") {$st = "u";}
  else {die "ERROR:  Unrecognized protonation flag: $ARGV[$argn]\n";}
  $nm = $ARGV[$argn+1];
  if (defined $nameindex{$nm}) {
      $i = $nameindex{$nm};
      push (@endlist, $i);
      $endstate[$i] = $st;
      }
  else {die "No such site as $nm\n";}
}

for (++$argn; $argn<@ARGV && $ARGV[$argn] ne "-with" ; $argn += 2) {
  if ($ARGV[$argn] eq "-p") {$st = "p";}
  elsif ($ARGV[$argn] eq "-u") {$st = "u";}
  else {die "ERROR:  Unrecognized protonation flag: $ARGV[$argn]\n";}
  $nm = $ARGV[$argn+1];
  if (defined $nameindex{$nm}) {
      $i = $nameindex{$nm};
      push (@reflist, $i);
      $curstate[$i] = $st;
  }
  else {die "No such site as $nm\n";}
}

# NEED A CHECK TO MAKE SURE @reflist AND @endlist ARE THE SAME
if (@endlist != @reflist) {die "endlist and reflist have different lengths\n";}
@sorted_endlist = sort @endlist;
@sorted_reflist = sort @reflist;
for ($j=0; $j<@sorted_endlist; ++$j) {
  if ($sorted_endlist[$j] != $sorted_reflist[$j]) {$notsame = 1;}
}
if ($notsame) {
  print STDERR "endlist and reflist don't have the same sites\n";
  print STDERR "     endlist:";
  foreach $j (@sorted_endlist) {print STDERR " $sitename[$j]";}
  print STDERR "\n     reflist:";
  foreach $j (@sorted_reflist) {print STDERR " $sitename[$j]";}
  die "\n";
}

# The remaining args should be the fixlist

for (++$argn; $argn<@ARGV; $argn+=2) {
  if ($ARGV[$argn] eq "-p") {$st = "p";}
  elsif ($ARGV[$argn] eq "-u") {$st = "u";}
  else {die "ERROR:  Unrecognized protonation flag: $ARGV[$argn]\n";}
  @adds = 0;
  @adds = grep (/^$ARGV[$argn+1]$/, @sitename);
  @adds || print "WARNING: No known residues match ${ARGV[$argn+1]}\n";
  if (grep ($_ eq $ARGV[$argn+1], @groups)) {
    warn "WARNING:  You have specified ${ARGV[$argn+1]} as fixed twice.\n";
  }
  else {
    push (@groups, $ARGV[$argn+1]);
  }
  ADDLOOP:
  foreach $nm (@adds) {
    $idx = $nameindex{$nm};
    # Don't let a reference state site get into the fixed list
    if (grep ($idx == $_, @reflist)) {
      print "NOTE: reference site $nm excluded from fixed list\n";
      next ADDLOOP;
    }
    # Redundecy check
    if (grep($idx == $_, @fixed)) {
      print "NOTE: $nm specified more than once as fixed site,\n",
             "last spec takes precedence.\n";
    }
    else {
      push (@fixed, $idx);
    }
    if ($groupname[$idx] !~ /\b$ARGV[$argn+1]\b/) {
      $groupname[$idx] .= "$ARGV[$argn+1] ";
    }
    $curstate[$idx] = $st;
  }
}


$total  = 0;
$totalpHfac = 0;
foreach $iref (@reflist) {
  while ($curstate[$iref] ne $endstate[$iref]) {

    if ($curstate[$iref] eq "u" && $endstate[$iref] eq "p") {
      $direction = "p";
    }
    else {
      $direction = "u";
    }
  
    # Terms due to pH and intrinsic pK
    if ($direction eq "p") {
      $pkintdiff = -$pkint[$iref];
      $pHfac = 1;
      print "Protonating site ${sitename[$iref]}:\n";
      print "     Intrisic: pkintdiff = $pkintdiff, pH term = $pHfac * pH\n";
    }
    else {
      $pkintdiff = $pkint[$iref];
      $pHfac = -1;
      print "Deprotonating site ${sitename[$iref]}:\n";
      print "     Intrisic: pkintdiff = $pkintdiff, pH term = $pHfac * pH\n";
    }
    $cumshift[$iref] = $pkintdiff;
  
    # Terms due to interaction with other sites
    @curfixed = (@reflist, @fixed);
    foreach $gr (@groups) {$groupshift{$gr} = 0;}
    CURFIXLOOP:
    foreach $ifix (@curfixed) {
      if ($ifix == $iref) {
	print "     Skipping same site ${sitename[$iref]}\n";
	next CURFIXLOOP;
      }
      # Pretend we're protonating site iref, then fix up sign later ..
      # Protonating a cationic site adds positive charge (raises protonation e)
      $shft = 0;
      if ($curstate[$ifix] eq "p" && $c_or_a[$ifix] eq "C") {
	$shft = $g[$iref*$numsites + $ifix]*$conv;
      }
      # Deprotonating an anionic site adds a negative charge (lowers prot. e)
      elsif ($curstate[$ifix] eq "u" && $c_or_a[$ifix] eq "A") {
	 $shft = -$g[$iref*$numsites + $ifix]*$conv;
      }
      # Other cases (prot. anionic, deprot. cationic) are already
      # in the "standard background" so no shift is calculated
  
      # Now fix up the sign
      if ($direction ne "p") {$shft = -$shft;}
      @ingroups = split(' ', $groupname[$ifix]);
      foreach $gr (@ingroups) {
	$groupshift{$gr} += $shft;
      }
      if (grep($ifix==$_, @reflist)) { # Print its contr if its in reflist ..
	print "     Shift due to ref site $sitename[$ifix] = $shft\n";
      }
      $cumshift[$iref] += $shft;
    }
    foreach $gr (@groups) {
      print "     Shift due to site group $gr = ${groupshift{$gr}}\n";
    }
    print "Free energy of this step = ${cumshift[$iref]} + pH * ($pHfac)\n\n";
    $curstate[$iref] = $direction;
    $total += $cumshift[$iref];
    $totalpHfac += $pHfac;
  }
}

print "TOTAL Free energy change = $total + pH * ($totalpHfac)\n";

