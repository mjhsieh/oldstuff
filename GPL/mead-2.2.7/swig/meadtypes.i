/* meadtypes.i --- common types used by MEAD -*- C++ -*-
 * ID: $Id: meadtypes.i,v 1.20 2002/02/03 18:23:56 bashford Exp $
 *
 * Author: John Bergsma <johnbergsma@hotmail.com>
 * Created: 2001/05/17 17:13:21
 */

%{
#include <strstream>
#include "stringobject.h"
%}

%{
#define SWIG_module_name "MEAD"
static PyObject *ThisModuleName;
%}

%init %{
    ThisModuleName=PyImport_ImportModule(SWIG_module_name);
    import_array();
%}

// Typemaps for strings
//
#define STRING_IN_MAP(who) %typemap who string * \
{ $target = new string (PyString_AsString ($source)); }

#define STRING_OUT_MAP(who) %typemap who string * \
{ $target = PyString_FromString ($source->c_str()); }

STRING_IN_MAP ((memberin))
STRING_IN_MAP ((python, in))
STRING_IN_MAP ((python, varin))

STRING_OUT_MAP ((memberout))
STRING_OUT_MAP ((python, out))
STRING_OUT_MAP ((python, varout))

%typemap(python,in) const string& ($basetype str) {
  if (PyString_Check($source)) {
    str = string(PyString_AsString($source));
    $target = &str;
  }
  else {
    PyErr_SetString(PyExc_TypeError, "arg must be a string");
    return NULL;
  }
}

%typemap(python,in) string ($basetype str) {
  if (PyString_Check($source)) {
    str = string(PyString_AsString($source));
    $target = &str;
  }
  else {
    PyErr_SetString(PyExc_TypeError, "arg must be a string");
    return NULL;
  }
}

/*
 * Use this macro when you need to make a pointer
 * to a Python shadow object (ie. FooPtr)
 * You just say:
 *
 * shadowClassPtr(Foo);
 *
 * in your interface file, and SWIG will generate a FooPtr() C function
 * for you. You can wrap your pointer like this:
 * PyObject *o=FooPtr((void *)source, 1); # 1 for own, 0 for the other case
 *
 * The argout typemaps can thus be written as
 *
 * %typemap(python,argout) Foo **{
 *    $target=t_output_helper($target,FooPtr((void*)*$source,1));
 * }
 * %typemap(python,ignore) Foo **(Foo *temp){
 *    $target=&temp;
 * }
 *
 */

%define shadowClassPtr(T)
%{
static PyObject *T##Ptr_o=NULL;
PyObject *T##Ptr(void *source,int own) {
    char cbuf[512];
    PyObject *o;
    if(T##Ptr_o==NULL)T##Ptr_o=PyObject_GetAttrString(ThisModuleName,"T""Ptr");
    SWIG_MakePtr(cbuf,(void*)(source),SWIGTYPE_p_##T);
    o=PyObject_CallFunction(T##Ptr_o, "(s)",cbuf);
    if(own){
       PyObject *val = PyInt_FromLong(1L);
       PyObject_SetAttrString(o, "thisown", val);
    }
    return o;
}
%}
%enddef

// Use these macros for stl typemaps

//
// Vector
//

// Typemaps for converting Python vector types to/from C++ vector<T> types
// To apply these typemaps define a typedef vector<T> as vector_T

#define Pystring_Check		PyString_Check
#define Pystring_Asstring	PyString_AsString
#define Pyint_Check		PyInt_Check
#define Pyint_Asint		PyInt_AsLong
#define Pylong_Check		PyLong_Check
#define Pylong_Aslong		PyInt_AsLong
#define Pyfloat_Check		PyFloat_Check
#define Pyfloat_Asfloat		PyFloat_AsDouble
#define Pydouble_Check		PyFloat_Check
#define Pydouble_Asdouble	PyFloat_AsDouble

// This is for basic Python types that have explicit conversion functions:
//    string
//    int
//    long
//    float
//    double
//
#define VECTOR_BASICTYPE_IN(T) %typemap(python,in) vector_##T & ($basetype vpe) { \
  if (PyList_Check($source)) { \
    int size = PyList_Size($source); \
    for (int i=0; i<size; ++i) { \
      PyObject *o = PyList_GetItem($source, i); \
      if (Py##T##_Check(o)) { \
        vpe.push_back( T(Py##T##_As##T (o))); \
      } \
      else { \
        PyErr_SetString(PyExc_TypeError, "Conversion of vector element to T failed"); \
        return NULL; \
      } \
    } \
    $target = &vpe; \
  } \
  else { \
    PyErr_SetString(PyExc_TypeError, "not a list"); \
    return NULL; \
  } \
}

// This is for wrapped types

#define VECTOR_IN(T) %typemap(python,in) vector_##T & ($basetype vpe) { \
  if (PyList_Check($source)) { \
    int size = PyList_Size($source); \
    for (int i=0; i<size; ++i) { \
      PyObject *o = PyList_GetItem($source, i); \
      T *peptr; \
      if ((SWIG_ConvertPtr(o,(void **) & peptr, SWIGTYPE_p_##T,1)) != -1) { \
        vpe.push_back(*peptr); \
      } \
      else { \
        PyErr_SetString(PyExc_TypeError, "Conversion of vector element to T failed"); \
        return NULL; \
      } \
    } \
    $target = &vpe; \
  } \
  else { \
    PyErr_SetString(PyExc_TypeError, "not a list"); \
    return NULL; \
  } \
}

#define VECTOR_PTR_IN(T) %typemap(python,in) vector_##T * ($basetype vpe) { \
  if (PyList_Check($source)) { \
    int size = PyList_Size($source); \
    for (int i=0; i<size; ++i) { \
      PyObject *o = PyList_GetItem($source, i); \
      T *peptr; \
      if ((SWIG_ConvertPtr(o,(void **) & peptr, SWIGTYPE_p_##T,1)) != -1) { \
        vpe.push_back(*peptr); \
      } \
      else { \
        PyErr_SetString(PyExc_TypeError, "Conversion of vector element to T failed"); \
        return NULL; \
      } \
    } \
    $target = &vpe; \
  } \
  else { \
    PyErr_SetString(PyExc_TypeError, "not a list"); \
    return NULL; \
  } \
}

#define VECTOR_ARGOUT(T) %typemap(python,argout) vector_##T * { \
  PyObject *newlist = PyList_New(0); \
  if (!newlist) return NULL; \
  for (vector_##T##::const_iterator p = $source->begin(); p!=$source->end(); ++p) { \
    T *fp = new T##(*p); \
    PyObject *optr = T##Ptr((void *) fp, 1); \
    PyList_Append(newlist, optr); \
    Py_DECREF(optr); \
  } \
  PyList_SetSlice($arg, 0, PyList_Size($arg), newlist); \
  Py_DECREF(newlist); \
}

//
// List
//

// Typemaps for converting Python list types to/from C++ list<T> types
// To apply these typemaps define a typedef list<T> as list_T
//
// Note: These only convert wrapped types
//       If you want to convert basic types follow the example above for vectors

// Note this works on const list_T types
#define LIST_IN(T) %typemap(python,in) const list_##T & ($basetype vi) { \
  if (PyList_Check($source)) { \
    int size = PyList_Size($source); \
    for (int i=0; i<size; ++i) { \
      PyObject *o = PyList_GetItem($source, i); \
      T *a=0; \
      if ((SWIG_ConvertPtr(o,(void **) &a,SWIGTYPE_p_##T,1)) == -1) { \
        PyErr_SetString(PyExc_TypeError, \
                        "Conversion of list item to T failed"); \
        return NULL; \
      } \
      vi.push_back(*a); \
    } \
    $target = &vi; \
  } \
  else { \
    PyErr_SetString(PyExc_TypeError, "not a list"); \
    return NULL; \
  } \
}

// To apply define the argout name as OutValue
#define LIST_ARGOUT(T) %typemap(python,argout) list_##T##& OutValue { \
  PyObject *newlist = PyList_New(0); \
  if (!newlist) return NULL; \
  for (list_##T##::const_iterator p = $source->begin(); p!=$source->end(); ++p) { \
    T *fp = new T##(*p); \
    PyObject *optr = T##Ptr((void *) fp, 1); \
    PyList_Append(newlist, optr); \
    Py_DECREF(optr); \
  } \
  PyList_SetSlice($arg, 0, PyList_Size($arg), newlist); \
  Py_DECREF(newlist); \
}

// Use this if you want pointers to list elements passed back.
// T needs to a pointer type
#define LIST_PTR_OUT(T) %typemap(python,out) list_##T##Ptr * { \
  $target = PyList_New(0); \
  if (!$target) return NULL; \
  for (list_##T##Ptr::const_iterator p = $source->begin(); p!=$source->end(); ++p) { \
    PyObject *optr = T##Ptr((void *) (*p), 0); \
    PyList_Append($target, optr); \
    Py_DECREF(optr); \
  } \
  delete $source; \
}

// Use this if you want copies of list elements passed back
#define LIST_OUT(T) %typemap(python,out) list_##T * { \
  $target = PyList_New(0); \
  if (!$target) return NULL; \
  for (list_##T##::const_iterator p = $source->begin(); p!=$source->end(); ++p) { \
    T *fp = new T##(*p); \
    PyObject *optr = T##Ptr((void *) fp, 1); \
    PyList_Append($target, optr); \
    Py_DECREF(optr); \
  } \
  delete $source; \
}

// Use this if you want copies of list elements passed back
// Leaks memory FIXME!
#define LIST_OF_PAIR_OUT(T1, T2) %typemap(python,out) list_##T1##_##T2##_pair * { \
  $target = PyList_New(0); \
  if (!$target) return NULL; \
  for (list_##T1##_##T2##_pair::const_iterator p = $source->begin(); p!=$source->end(); ++p) { \
    T1 *fp = new T1##(p->first); \
    PyObject *optr1 = T1##Ptr((void *) fp, 1); \
    T2 *sp = new T2##(p->second); \
    PyObject *optr2 = T2##Ptr((void *) sp, 1); \
    PyObject *newtuple = PyTuple_New(2); \
    PyTuple_SetItem(newtuple, 0, optr1); \
    PyTuple_SetItem(newtuple, 1, optr2); \
    PyList_Append($target, newtuple); \
  } \
  delete $source; \
}

// Use this if you want a list of pairs passed back where T2 is a pointer
// Leaks memory FIXME!
#define LIST_OF_PAIR_PTR_OUT(T1, T2) %typemap(python,out) list_##T1##_##T2##_pairPtr * { \
  $target = PyList_New(0); \
  if (!$target) return NULL; \
  for (list_##T1##_##T2##_pairPtr::const_iterator p = $source->begin(); p!=$source->end(); ++p) { \
    T1 *fp = new T1##(p->first); \
    PyObject *optr1 = T1##Ptr((void *) fp, 1); \
    PyObject *optr2 = T2##Ptr((void *) &(*p->second), 0); \
    PyObject *newtuple = PyTuple_New(2); \
    PyTuple_SetItem(newtuple, 0, optr1); \
    PyTuple_SetItem(newtuple, 1, optr2); \
    PyList_Append($target, newtuple); \
  } \
  delete $source; \
}


// Map

// Typemaps for converting Python dict types to/from C++ map<T1, T2> types
// To apply these typemaps define a typedef map<T1,T2> as map_T1_T2
//
// Note: These only convert wrapped types
//       If you want to convert basic types follow the example above for vectors
//
#define MAP_IN(T1, T2) %typemap(python,in) map_##T1##_##T2 & ($basetype vi) { \
  if (PyDict_Check($source)) { \
    PyObject *key, *value; \
    int pos = 0; \
    while (PyDict_Next($source, &pos, &key, &value)) { \
      T1 *a=0; \
      if ((SWIG_ConvertPtr(key,(void **) &a,SWIGTYPE_p_##T1,1)) == -1) { \
        PyErr_SetString(PyExc_TypeError, \
                        "Conversion of map key to T1 failed"); \
        return NULL; \
      } \
      T2 *b=0; \
      if ((SWIG_ConvertPtr(value,(void **) &a,SWIGTYPE_p_##T2,1)) == -1) { \
        PyErr_SetString(PyExc_TypeError, \
                        "Conversion of map value to T2 failed"); \
        return NULL; \
      } \
      vi.insert(map_##T1##_T2##::value_type(*a,*b)); \
    } \
    $target = &vi; \
  } \
  else { \
    PyErr_SetString(PyExc_TypeError, "not a dict"); \
    return NULL; \
  } \
}

// To apply define the argout name as OutValue
#define MAP_ARGOUT(T1, T2) %typemap(python,argout) map_##T1##_##T2##& OutValue { \
  PyDict_Clear($arg); \
  for (map_##T1##_##T2##::const_iterator p = $source->begin(); p!=$source->end(); ++p) { \
    T1 *key = new T1##(p->first); \
    PyObject *optr1 = T1##Ptr((void *) key, 1); \
    T2 *value = new T2##(p->second); \
    PyObject *optr2 = T2##Ptr((void *) value, 1); \
    PyDict_SetItem($arg, optr1, optr2); \
    Py_DECREF(optr1); \
    Py_DECREF(optr2); \
  } \
}

%{
#include "Numeric/arrayobject.h"
%}

// Map to a 3-D Python Numeric Array.
// To use call a function that returns a Numeric_Array_T * generated from
// a CubeLatSpec&. Remember to convert the array indicies to Fortran style.
#define NUMERIC_3D_ARRAY_OUT(T) %typemap(python,out) Numeric_3D_Array_##T * { \
  if ($source != 0) { \
    int dims[3]; \
    dims[0] = dims[1] = dims[2] = arg1->get_grid_dim(); \
    PyArrayObject * array_data = (PyArrayObject *) PyArray_FromDimsAndData (3, dims, PyArray_##T, (char *) $source); \
    if (!array_data) return NULL; \
    array_data->flags |= OWN_DATA; \
    $target = (PyObject *) array_data; \
  } \
  else { \
    PyErr_SetString (PyExc_RuntimeError, "Numeric_Array_""T"": memory alloc failed"); \
    $target = 0; \
  } \
}

// Map to a 1-D Python Numeric Array from an arbitrary set of points
// To use call a function that returns a Numeric_Array_1D * from npts and data
#define NUMERIC_1D_ARRAY_OUT(T) %typemap(python,out) Numeric_1D_Array_##T * { \
  if ($source != 0) { \
    int dim = arg1; \
    PyArrayObject * array_data = (PyArrayObject *) PyArray_FromDimsAndData (1, &dim, PyArray_##T, (char *) $source); \
    if (!array_data) return NULL; \
    array_data->flags |= OWN_DATA; \
    $target = (PyObject *) array_data; \
  } \
  else { \
    PyErr_SetString (PyExc_RuntimeError, "Numeric_Array_""T"": memory alloc failed"); \
    $target = 0; \
  } \
}

// Input conversions for vector types
%{
#include <vector>
#include "MEAD/Coord.h"
typedef std::vector<Coord> vector_Coord;
typedef std::vector<string> vector_string;
typedef std::vector<int> vector_int;
typedef std::vector<float> vector_float;
typedef std::vector<double> vector_double;
%}
VECTOR_IN(Coord)
VECTOR_BASICTYPE_IN(string)
VECTOR_BASICTYPE_IN(int)
VECTOR_BASICTYPE_IN(float)
VECTOR_BASICTYPE_IN(double)

// Input conversion for list<Atom> == list_Atom
%{
#include <list>
#include "MEAD/Atom.h"
%}
shadowClassPtr(Atom)
LIST_IN(Atom)

// Input conversion for map<AtomID,Atom> == map_AtomID_Atom
%{
#include <map>
#include "MEAD/AtomID.h"
typedef std::map<AtomID,Atom> map_AtomID_Atom;
%}
shadowClassPtr(AtomID)
MAP_IN(AtomID,Atom)

// For constructing an output list of map::value_type pairs
// Leaks memory FIXME!
//%{
//typedef std::pair<AtomID,Atom> AtomID_Atom_pair;
//typedef std::list<AtomID_Atom_pair> list_AtomID_Atom_pair;
//%}
//LIST_OF_PAIR_OUT(AtomID,Atom)

// For constructing an output list of AtomSet map::keys
%{
typedef std::list<AtomID> list_AtomID;
%}
LIST_OUT(AtomID)

// For constructing an output list of pointers to AtomSet map::values
%{
typedef std::list<Atom*> list_AtomPtr;
%}
LIST_PTR_OUT(Atom)

// In/Out conversions for list<PointCharge> = list_PointCharge
%{
#include "MEAD/PointCharge.h"
%}
shadowClassPtr(PointCharge)
LIST_IN(PointCharge)
LIST_OUT(PointCharge)

// Output conversion for list<CubeLatSpec> = list_CubeLatSpec
%{
#include "MEAD/CubeLatSpec.h"
typedef std::list<CubeLatSpec> list_CubeLatSpec;
%}
shadowClassPtr(CubeLatSpec)
LIST_OUT(CubeLatSpec)

// Output conversion for vector<Legendre> = vector_Legendre
%{
#include "MEAD/SphericalHarmonic.h"
typedef std::vector<Legendre> vector_Legendre;
%}
shadowClassPtr(Legendre)
VECTOR_PTR_IN(Legendre)
VECTOR_ARGOUT(Legendre)

// Output conversions to Python Numeric Arrays
%{
typedef float Numeric_3D_Array_FLOAT;
typedef int Numeric_3D_Array_INT;
typedef float Numeric_1D_Array_FLOAT;
typedef int Numeric_1D_Array_INT;
%}
NUMERIC_3D_ARRAY_OUT(INT)
NUMERIC_3D_ARRAY_OUT(FLOAT)
NUMERIC_1D_ARRAY_OUT(INT)
NUMERIC_1D_ARRAY_OUT(FLOAT)

// Don's hand-rolled typemap for x by 4 numeric arrays
%typemap(python,in) Numeric_Nby4_DOUBLE * {
  if (PyArray_Check($source)) {
    PyArrayObject *array = (PyArrayObject *) $source;
    if (array->nd != 2) {
      PyErr_SetString(PyExc_TypeError,
		      "array must be two-diminsional");
      return NULL;
    }
    if (array->dimensions[1] != 4) {
      PyErr_SetString(PyExc_TypeError,
		      "array's second diminsion size must be 4");
      return NULL;
    }
    if (array->descr->type_num != PyArray_DOUBLE) {
      PyErr_SetString(PyExc_TypeError,
		      "array contents must be of type double");
      return NULL;
    }
    $target = array;
  }
  else {
    PyErr_SetString(PyExc_TypeError, "not an array");
    return NULL;
  }
}
