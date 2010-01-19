#ifndef _ioNDIM_h
#define _ioNDIM_h

/*

Useful general-purpose i/o operators go here

*/

#include <iostream>
#include <fstream>
#include "kelp.h"
#include "GridNDIM.h"

#undef OUTPUT
#define OUTPUT(X) { if (mpLeader()) cout << X << flush; }

ostream& operator << (ostream& o, const GridNDIM<double> &G);

#endif
