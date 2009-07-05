#!/usr/bin/perl -w
# Title:  A script to retreive information from installation records
#         in Mac OS X Leopard or greater.
# Author: Mengjuei Hsieh
# 
#
use strict;
use warnings;
use Foundation;

my $receipts_dir = '/Library/Receipts';
my %package_hash = ();
my $bool_found = 0;

sub main {
    if ( @ARGV < 1 ) {
        rtfm();
    } else {
        get_all_package_ids();
    }
    if ( $ARGV[0] eq "list" ) {
        foreach my $package_id (sort keys %package_hash) {
            if ( $package_hash{$package_id} eq "bom" ) {
               print "old-styled $package_id\n";
            } elsif ( $package_hash{$package_id} eq "flat" ) {
               print "new-styled $package_id\n";
            }
        }
    } elsif ( ($ARGV[0] eq "info") and (@ARGV > 1) ) {
        my $package_id = ();
        foreach $package_id (keys %package_hash) {
            if ( $ARGV[1] eq $package_id ) {$bool_found = 1;}
        }
        if ($bool_found == 1) {
            get_package_info($ARGV[1]);
        } else {
            die "No such package: $ARGV[1]\n";
        }
    } elsif ( ($ARGV[0] eq "files") and (@ARGV > 1) ) {
        my $package_id = ();
        foreach $package_id (keys %package_hash) {
            if ( $ARGV[1] eq $package_id ) {$bool_found = 1;}
        }
        if ($bool_found == 1) {
            get_package_files($ARGV[1]);
        } else {
            die "No such package: $ARGV[1]\n";
        }
    } else {
        rtfm();
    }
}

sub get_package_files {
    my $pkg_id = shift;
    if ( $bool_found == 0 ) {
        die "incorrect use of get_package_files\n";
    } elsif ( $package_hash{$pkg_id} eq "bom" ) {
        my $package_path = "$receipts_dir/$pkg_id";
        my $package_bom_path = undef;
        my $package_name = "$pkg_id";
        $package_name =~ s/\.pkg$//;
        foreach my $package_bom_candidate (
            "$package_path/Contents/Archive.bom",
            "$package_path/Contents/Resources/Archive.bom",
            "$package_path/Contents/Resources/$package_name.bom",
        ) {
            if ((-f $package_bom_candidate) && (-r $package_bom_candidate)) {
                $package_bom_path = $package_bom_candidate;
                last;
            }
        }
        foreach my $pkg_files (execute_command(['lsbom'],['-f','-l','-d','-p','f',$package_bom_path])){
            $pkg_files =~ s/^\.\///;
            print "$pkg_files\n";
        }
    } elsif ( $package_hash{$pkg_id} eq "flat" ) {
        foreach my $pkg_files (execute_command(['pkgutil'],['--files',$ARGV[1]])) {
            print "$pkg_files\n";
        }
    }
}

sub get_package_info {
    my $pkg_id = shift;
    my $myobj = ();
    if ( $bool_found == 0 ) {
        die "incorrect use of get_package_info\n";
    } elsif ( $package_hash{$pkg_id} eq "bom" ) {
        my $Infofile = $receipts_dir . "/" . $pkg_id . "/Contents/Info.plist";
        my $bomplist = NSDictionary->dictionaryWithContentsOfFile_( $Infofile );
        if ( ! $bomplist or ! $$bomplist ) {
            print "One of the following failed: $bomplist - $$bomplist\n";
            exit 1;
        }
        $myobj = $bomplist->objectForKey_("CFBundleIdentifier");
        print "package-id: " . $myobj->description()->UTF8String() . "\n";
        $myobj = $bomplist->objectForKey_("CFBundleShortVersionString");
        print "version: " . $myobj->description()->UTF8String() . "\n";
        $myobj = $bomplist->objectForKey_("IFPkgFlagDefaultLocation");
        print "volume: " . $myobj->description()->UTF8String() . "\n";
        $myobj = $bomplist->objectForKey_("IFPkgRelocatedPath");
        print "location: " . $myobj->description()->UTF8String() . "\n";
    } elsif ( $package_hash{$pkg_id} eq "flat" ) {
        foreach my $package_info (execute_command(['pkgutil'],['--pkg-info',$ARGV[1]])) {
            print "$package_info\n";
        }
    }
}

sub rtfm {
    prn_usage();
    exit 0;
}

sub get_all_package_ids {
    if (opendir(RECEIPTS_DIR,"$receipts_dir")) {
        while (defined(my $package_id = readdir(RECEIPTS_DIR))) {
            if ($package_id =~ m|^.*\.pkg$|) {
                 $package_hash{$package_id} = 'bom';
            }
        }
        closedir(RECEIPTS_DIR);
    }
    foreach my $package_id (execute_command(['pkgutil'],['--pkgs'])) {
        $package_hash{$package_id} = 'flat';
    }
}

sub prn_usage {
    print "Usage: ${0} <verb> [option]\n";
    print "     <verb> is one of the following:\n";
    print "     list                  (List all currently installed package with records.)\n";
    print "     info PKGID            (Show metadata about PKGID)\n"; 
    print "     files PKGID           (List files installed by the specified package)\n"; 
}

sub execute_command {
    my $command_prog_ref = shift;
    my $command_argv_ref = shift;
    my @result_output = ();
    if (open(CHILD, "-|")) {
        while (<CHILD>) {
            chomp;
            push(@result_output,$_);
        }
        close(CHILD);
    } else {
        exec(@$command_prog_ref, @$command_argv_ref);
    }

    return @result_output;
}

main();
