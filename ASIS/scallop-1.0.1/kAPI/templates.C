#include <stdio.h>

#if defined(__GNUC__)
/* Templates are explicitly instantiated for g++.  
   The other options scare me. */

#include "msg.h"
#include "msg.C"
#include "List.h"
#include "List.C"
#include "ArrayNDIM.h"
#include "ArrayNDIM.C"
#include "GridNDIM.h"
#include "GridNDIM.C"
#include "XArrayNDIM.h"
#include "XArrayNDIM.C"
#include "s_gridNDIM.h"
#include "VectorMoverNDIM.h"
#include "VectorMoverNDIM.C"

template mpMessageID * mpARecv<char>(int, char *, int);
template mpMessageID * mpASend<double>(int, double const *, int, int);
template class KeLP::ListLink<MotionItemNDIM *>;
template class KeLP::ListIterator<MotionItemNDIM *>;
template class KeLP::List<MotionItemNDIM *>;
template class ArrayNDIM<double>;
template class VectorMoverNDIM<s_gridNDIM, double>;
template class XArrayNDIM<s_gridNDIM>;
template class GridNDIM<double>;

#endif // defined(__GNUC__)

#if defined(__xlC__)

#include "s_gridNDIM.h"
#include "XArrayNDIM.h"
#include "VectorMoverNDIM.h"

#pragma define(XArrayNDIM<s_gridNDIM>)
#pragma define(VectorMoverNDIM<s_gridNDIM,double>)

#include "XArrayNDIM.C"
#include "VectorMoverNDIM.C"

#endif // defined(__xlC__)
