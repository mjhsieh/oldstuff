#!/bin/sh -x

molname=$1

if [ -f $molname.pqr -a -f $molname.ogm -a -f $molname.mgm ] ; then
  for taut in del eps del1 del2 del3 del4 del5 del6 del7 del8 del9 del10 del11
  do
    if [ -f $molname.sites.$taut ]
    then
      if [ -f $molname.sites ] ; then
        /bin/rm -f $molname.sites
      fi
      ln -s $molname.sites.$taut $molname.sites
      time ./multiflex -epsave_oldway -epsin 4.0 -blab1 -ionicstr 0.1 \
           $molname > $molname-$taut.el.out
      mv $molname.pkint $molname-$taut.pkint
      mv $molname.g $molname-$taut.g
    else 
      echo "ERROR $molname.sites.$taut non found.  Continuing..."
    fi
  done
else
  echo "ERROR: molname, $molname, is wrong or files are missing"
fi

