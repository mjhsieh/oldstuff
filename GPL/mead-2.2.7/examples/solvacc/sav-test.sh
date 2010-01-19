#!/bin/sh
BINDIR=${BINDIR:=..}
set echo
${BINDIR}/sav-driver -blab2 trna
${BINDIR}/sav-driver -dimension 41 -spacing 0.5 trna
${BINDIR}/sav-driver -dimension 41 -spacing 0.1 trna
${BINDIR}/sav-driver -dimension 61 -spacing 1.0 trna
${BINDIR}/sav-driver -dimension 61 -spacing 0.5 trna
${BINDIR}/sav-driver -dimension 61 -spacing 0.1 trna
${BINDIR}/sav-driver -dimension 81 -spacing 1.0 trna
${BINDIR}/sav-driver -dimension 81 -spacing 0.5 trna
${BINDIR}/sav-driver -dimension 81 -spacing 0.1 trna
${BINDIR}/sav-driver -dimension 101 -spacing 1.0 trna
${BINDIR}/sav-driver -dimension 101 -spacing 0.5 trna
${BINDIR}/sav-driver -dimension 101 -spacing 0.1 trna
${BINDIR}/sav-driver -dimension 121 -spacing 1.0 trna
${BINDIR}/sav-driver -dimension 121 -spacing 0.5 trna
${BINDIR}/sav-driver -dimension 121 -spacing 0.1 trna
unset echo
