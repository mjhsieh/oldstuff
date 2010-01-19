#!/bin/sh -vx
molname=COx-AmberBondi
  ./run_mol_multimead.sh $molname
  ./make-globals.pl  $molname
  ./runmcti.sh $molname
  ./collect-curves.pl mcti.out.global > curves.out
  grep 'pK(1/2)' curves.out > pkhalf.out
