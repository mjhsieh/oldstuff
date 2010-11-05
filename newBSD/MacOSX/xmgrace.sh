myMAXPATH=50000
if [ "${1}test" != "test" ]; then
   myFILEPATH="$(/usr/bin/dirname "$1")"
   if [ "${myFILEPATH}test" != "test" ]; then
      cd "$myFILEPATH"
   fi
   /opt/local/bin/xmgrace -g 1130x894 -maxpath $myMAXPATH "$1"
else
   /opt/local/bin/xmgrace -g 1130x894 -maxpath $myMAXPATH
fi
