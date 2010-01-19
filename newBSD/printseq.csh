#!/bin/csh -f
# Only single chain
switch ( $#argv )
  case 1:
    set tmpseq=`cat $1 | grep "^ATOM" | cut -c 18-20,21-26 | uniq | awk '{print $1,$2}'`
    breaksw
  case 2:
    set range=`echo $2 | sed -e 's/-/ /'`
    setenv range1 $range[1]
    setenv range2 $range[2]
    set tmpseq=`cat $1 | grep "^ATOM" | cut -c 18-20,21-26 | uniq \
	| awk '{if ($2 >= ENVIRON["range1"] && $2 <= ENVIRON["range2"] ) print $1,$2}'`
    breaksw
  default:
    echo error N
    exit 1
    breaksw
endsw

echo $tmpseq | sed -e 's/[:A-Z:]/ /g' | awk '{print $1,$NF}'
echo $tmpseq \
	| sed -e 's/[:0-9:]//g;\
		s/ALA/A/g;\
		s/CYS/C/g;s/CYX/C/g;\
		s/ASP/D/g;\
		s/GLU/E/g;\
		s/PHE/F/g;\
		s/GLY/G/g;\
		s/HIS/H/g;s/HID/H/g;s/HSD/H/g;\
		s/ILE/I/g;\
		s/LYS/K/g;\
		s/LEU/L/g;\
		s/MET/M/g;\
		s/ASN/N/g;\
		s/PRO/P/g;\
		s/GLN/Q/g;\
		s/ARG/R/g;\
		s/SER/S/g;\
		s/THR/T/g;\
		s/TRP/W/g;\
		s/VAL/V/g;\
		s/TYR/Y/g;' \
	| sed -e 's/ //g;' \
	| tr "[:upper:]" "[:lower:]"
