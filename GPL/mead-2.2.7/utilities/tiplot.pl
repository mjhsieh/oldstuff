#!/usr/bin/perl

# This is a replacement for curdisp.  It derives from the tiplot
# that I wrote during the mb project to process
# Dave's "pkhalf.out" files, which was derived from my original sh/awk
# program curdisp.
# It's purpose is to display and optionally plot titration curves.
# It starts up a gnuplot and
# displays titration curves for each site in turn, asking if you want
# to print it.  Leaves behind a "curves.ps" or something you can send
# to printer for hardcopy.

# Usage: tiplot molname

# Where molname will appear in the titles.

$prefix = $ARGV[0];
if (! $prefix) {
  die "tiplot: ERROR: no molname prefix given on command line.\n";
  }

open (PKHALF, "${prefix}pkout")
    || open (PKHALF, "${prefix}.pkout")
    ||  die "Couldn't open ${prefix}pkout\n";

# Open up two gnuplot pipes -- one for interactively showing all
# curves and another for plotting selected ones.
open (SHOWPLOT, "| gnuplot")  || die "Couldn't start a gnuplot pipe\n";
open (PRINPLOT, "| gnuplot")  || die "Couldn't start a gnuplot pipe\n";

# Make the gluplot pipes "hot" ones.
# YOW! see p 110 of the perl book for these one-liners!
select((select(SHOWPLOT), $| = 1)[0]);
select((select(PRINPLOT), $| = 1)[0]);

print SHOWPLOT "set term x11\n";

print PRINPLOT "set term postscript\n";
$psplotfile = "curves.ps";
print PRINPLOT "set output \"$psplotfile\"\n";

$sitename = "";
$datfile= "";
$ptcount = 0;
$samples = 160;    # The gnuplot default sample rate
$num_printed = 0;

line: while (<PKHALF>) {
#    chop;	# strip record separator
    @Fld = split(' ');
    if (/^Site/ || /^Whole protein/) {
        # Plot the previous site if need be.
        if ($sitename) {
	    print "Here with site $sitename\n";
            &PlotIt;
            }
        # Prepare for a new site.
        $sitename = $Fld[1];
        $datfile = $sitename . ".dat";
        open (DAT, ">$datfile") || die "Couldn't open data output file, ",
                                       "$datfile\n";
        $ptcount = 0;
        next line;
        }
    # Process a possible data point
    if ($Fld[0] =~ /^-?[0-9]/ && $Fld[1] =~ /^-?[0-9]/ && $sitename) {
        print DAT $_;
	++$ptcount;
    }
    elsif ($sitename) {
        &PlotIt;
        $sitename = "";
        $ptcount = 0;
        }

    }

# Take care of the last one...
  & PlotIt;

if ($num_printed > 0) {
    print "The $num_printed curves requested for printing are in ".
	"Postscript file, $psplotfile\n";
}
else {
    unlink $psplotfile;
}

        # Close the previous sites curve file, and display it.
sub PlotIt {
            close DAT;
            if ($ptcount) {
                if ($ptcount > $samples) {
                    # Need to up the gnuplot sampling rate
                    $samples = $ptcount + 1;
                    print SHOWPLOT "set samples $samples\n";
                    print PRINPLOT "set samples $samples\n";
                    }
                print SHOWPLOT "plot \"$datfile\" ",
                               "title '$prefix Site $sitename' with lines\n";
                print "Do you want to print this one?  [y/n]";
                $yesno = <STDIN>;
                if ($yesno =~ /[Yy][Ee]?[Ss]?/) {
		    ++$num_printed;
                    print PRINPLOT "plot \"$datfile\" ",
                                "title '$prefix Site $sitename' with lines\n";
                    sleep 2;
                    # Above is a gross hack to stop the unlink below from
                    # happening before PRINPLOT (a gnuplot proc) is done
                    # with the file.
                    }
                }
            else {
                print "Site $sitename doesn't titrate.  Skipping to next\n";
                print SHOWPLOT "clear\n";
                }
            unlink $datfile;
}
