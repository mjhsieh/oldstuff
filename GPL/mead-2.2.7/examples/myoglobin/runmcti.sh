#!/bin/sh -vx

# Run Paul Beroza's mcti program.

time ./mcti << EOF > mcti.out.global
$1.global.pkint
$1.global.g
1000		! number of full MC steps
5000		! number of reduced MC steps
-4.0		! starting pH
15.0		! final pH
0.1		! pH increment
2.0		! min_wint for pairs
0.000001	! redmc toler
0		! r.n. seed (0 = auto)
1		! 0=full, 1=reduced m.c.
mcti.log.global
EOF

