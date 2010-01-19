#!/bin/sh
myMAXPATH=50000
if [ "${1}test" != "test" ]; then
   cd "$(/usr/bin/dirname "$1")"
   /opt/local/bin/xmgrace -g 1158x808 -maxpath ${myMAXPATH} "$1"
else
   /opt/local/bin/xmgrace -g 1158x808 -maxpath ${myMAXPATH}
fi
