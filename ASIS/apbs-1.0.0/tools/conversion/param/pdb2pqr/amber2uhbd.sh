#!/bin/bash

# This converts AMBER format parameter files to UHBD.  Take just the numbered
# stuff for each residue, i.e.:
#
#  1   DUMM  DU    M    0  -1  -2    0.00      0.00      0.00       0.0000
#  2   DUMM  DU    M    1   0  -1    1.00      0.00      0.00       0.0000
#  3   DUMM  DU    M    2   1   0    1.00     90.00      0.00       0.0000
#  4   H5T   HO    M    3   2   1    1.00    120.00    180.00       0.4422
#  5   O5'   OH    M    4   3   2    0.96    101.43    -98.89      -0.6318
#  6   C5'   CT    M    5   4   3    1.44    119.00    -39.22      -0.0069
#  7   H5'1  H1    E    6   5   4    1.09    109.50     60.00       0.0754
#  8   H5'2  H1    E    6   5   4    1.09    109.50    -60.00       0.0754
#  9   C4'   CT    M    6   5   4    1.52    110.00    180.00       0.1629
# 10   H4'   H1    E    9   6   5    1.09    109.50   -200.00       0.1176
# 11   O4'   OS    S    9   6   5    1.46    108.86    -86.31      -0.3691
# 12   C1'   CT    B   11   9   6    1.42    110.04    105.60       0.0431
# 13   H1'   H2    E   12  11   9    1.09    109.50   -240.00       0.1838
# 14   N9    N*    S   12  11   9    1.52    109.59   -127.70      -0.0268
# 15   C8    CK    B   14  12  11    1.37    131.20     81.59       0.1607
# 16   H8    H5    E   15  14  12    1.08    120.00      0.00       0.1877
# 17   N7    NB    S   15  14  12    1.30    113.93    177.00      -0.6175
# 18   C5    CB    S   17  15  14    1.39    104.00      0.00       0.0725
# 19   C6    CA    B   18  17  15    1.40    132.42    180.00       0.6897
# 20   N6    N2    B   19  18  17    1.34    123.50      0.00      -0.9123
# 21   H61   H     E   20  19  18    1.01    120.00    180.00       0.4167
# 22   H62   H     E   20  19  18    1.01    120.00      0.00       0.4167
# 23   N1    NC    S   19  18  17    1.34    117.43    180.00      -0.7624
# 24   C2    CQ    B   23  19  18    1.33    118.80      0.00       0.5716
# 25   H2    H5    E   24  23  19    1.08    120.00    180.00       0.0598
# 26   N3    NC    S   24  23  19    1.32    129.17      0.00      -0.7417
# 27   C4    CB    E   26  24  23    1.35    110.80      0.00       0.3800
# 28   C3'   CT    M    9   6   5    1.53    115.78   -329.11       0.0713
# 29   H3'   H1    E   28   9   6    1.09    109.50     30.00       0.0985
# 30   C2'   CT    B   28   9   6    1.53    102.80    -86.30      -0.0854
# 31   H2'1  HC    E   30  28   9    1.09    109.50    120.00       0.0718
# 32   H2'2  HC    E   30  28   9    1.09    109.50    240.00       0.0718
# 33   O3'   OS    M   28   9   6    1.42    116.52   -203.47      -0.5232
# 
# and put it in a file labeled RES.prn, where RES is the residue name.

for file in *.prn; do
    resi=${file%.*}
    cat $file | awk -v x=$resi '{print x, $2, $11, 0.000, 0.000, $8}'
done

