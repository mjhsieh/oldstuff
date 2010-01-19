#!/bin/sh

# ./configure --enable-wrapping
# make MEAD_shadow.cc

echo WARNING: This file has bit-rotted.
echo It needs to be updated before use.
echo '   --ttn, 2000/06/12 17:14:43'
echo Terminating failurefully.
exit 1


mwddir=`pwd`/../mead-windows-distro
mkdir $mwddir
mkdir $mwddir/libmead
mkdir $mwddir/swig
mkdir $mwddir/apps
mkdir $mwddir/apps/solvate
mkdir $mwddir/gnustuff

for d in libmead apps/solvate swig gnustuff
    do
	backdir=`pwd`
	cd $d
	for ccf in *.cc *.h
	    do
		cppf=${mwddir}/${d}/`echo $ccf | sed 's/\.cc$/.cpp/'`
		echo MAKING $cppf
		sed 's/<String.h>/<GString.h>/' < $ccf >$cppf
	    done
        cd $backdir
    done

cat <<EOF > $mwddir/gnustuff/strstream.h
#include <strstrea.h>

EOF

/bin/rm -f $mwddir/gnustuff/String.h
cp gnustuff/String.hin $mwddir/gnustuff/GString.h
cp swig/MEAD.py $mwddir/swig/MEAD.py
