#!/usr/bin/env perl
$numArgs = $#ARGV + 1;
if($numArgs != 2){ die "Usage: removedefs.pl file.py defs\n"; }
$pyfile=shift @ARGV;
$deffile=shift @ARGV;

open(DEFFILE,$deffile);
while(<DEFFILE>){
   if(/^def\s+(.*)\(\*args\):$/o){
      $defs{$1} = 0;
   }else{
      die "Error: $deffile contains non-def entry $1 at line: $. \n"; }
}
close(DEFFILE);

$founddef=0;
open(PYFILE,$pyfile);
while(<PYFILE>){
  if($founddef){ --$founddef; next; }
  foreach $d (keys %defs){
     if(!$defs{$d} && /^def\s+$d\(\*args\):$/){
        $founddef=3;
        ++$defs{$d};
        last;
     }
  }
  if ($founddef){ next; }
  print;
}
close(PYFILE);

#print "\nFound the following def's\n";
foreach $d (keys %defs){
#  print "$d $defs{$d}\n";
   if(!$defs{$d}){
      print STDERR "WARNING: $d not found in $pyfile \n";
   }
   if($defs{$d} > 1){
      print STDERR "WARNING: $d found $defs{$d} times in $pyfile \n";
   }
}
