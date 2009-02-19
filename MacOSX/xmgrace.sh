#!/bin/sh
if [ "${1}test" != "test" ]; then
   cd "$(/usr/bin/dirname "$1")"
fi
/opt/local/bin/xmgrace -g 1158x808 -maxpath 50000 "$1"
