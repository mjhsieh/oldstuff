#!/bin/bash
(cd amber11.cvs; cvs -z3 update -Pd)
if [ -d amber11.cvs -a -d amber11.svn ]; then
   rsync -avx --exclude 'bin' --exclude 'exe' --exclude 'CVS' --exclude '.svn' \
         --delete amber11.cvs/ amber11.svn/
fi
cd amber11.svn
for missedfile in $(svn st | grep -e '^!' | sed -e 's/!	//'); do
   if [ -d "${missedfile}" ] ; then
      svn rm "${missedfile}"
   fi
done
for missedfile in $(svn st | grep -e '^!' | sed -e 's/! //'); do
   svn rm "${missedfile}"
done
for extrafile in $(svn st | grep -e '^?' | sed -e 's/? //'); do
   if [ -d "${extrafile}" ] ; then
      svn add "${extrafile}"
   fi
done
for extrafile in $(svn st | grep -e '^?' | sed -e 's/? //'); do
   svn add "${extrafile}"
done
