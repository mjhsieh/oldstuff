#!/bin/sh
set echo
./sav-driver -blab2 trna
./sav-driver -dimension 41 -spacing 0.5 trna
./sav-driver -dimension 41 -spacing 0.1 trna
./sav-driver -dimension 61 -spacing 1.0 trna
./sav-driver -dimension 61 -spacing 0.5 trna
./sav-driver -dimension 61 -spacing 0.1 trna
./sav-driver -dimension 81 -spacing 1.0 trna
./sav-driver -dimension 81 -spacing 0.5 trna
./sav-driver -dimension 81 -spacing 0.1 trna
./sav-driver -dimension 101 -spacing 1.0 trna
./sav-driver -dimension 101 -spacing 0.5 trna
./sav-driver -dimension 101 -spacing 0.1 trna
./sav-driver -dimension 121 -spacing 1.0 trna
./sav-driver -dimension 121 -spacing 0.5 trna
./sav-driver -dimension 121 -spacing 0.1 trna
unset echo
