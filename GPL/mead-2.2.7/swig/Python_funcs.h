/* Template functions to emulate Python objects
   Copyright (C) 1990--1995 by Donald Bashford.

    This source code file is part of the MEAD (Macroscopic
    Electrostatics with Atomic Detail) package of objects and
    programs.  MEAD is free software; you can redistribute it and/or
    modify it under the terms of the GNU General Public License as
    published by the Free Software Foundation; either version 1, or
    (at your option) any later version.

    This program is distributed in the hope that it will be useful,
    but WITHOUT ANY WARRANTY; without even the implied warranty of
    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the GNU
    General Public License for more details.

    You should have received a copy of the GNU General Public License
    along with this program; see the file COPYING.  If not, write to
    the Free Software Foundation, Inc., 675 Mass Ave, Cambridge, MA
    02139, USA.

    Donald Bashford can be contacted by electronic mail by the address,
    bashford@scripps.edu, or by paper mail at Department of Molecular
    Biology, The Scripps Research Institute, 10666 North Torrey Pines
    Road, La Jolla, California 92037.

$Id: Python_funcs.h,v 1.3 2004/12/06 23:21:45 bashford Exp $
*/

#ifndef _Python_funcs_h
#define _Python_funcs_h 1

#include <list>
#include <map>

using std::list;
using std::map;
using std::pair;
using std::copy;
using std::back_inserter;

#ifdef __GNUG__
#pragma implementation
#endif

/*
 * Functions to emulate a Python list
 */

// Append an item to a list
template<class T> void list_append (list<T>* lptr, const T& item)
{
  lptr->push_back(item);
}
// Length of list
template<class T> int list___len__ (list<T>* lptr)
{
  return lptr->size();
}
// Add (copy) all items from list l to end of list
template<class T> void list_extend (list<T>* lptr, const list<T>& l)
{
  copy(l.begin(), l.end(), back_inserter(*lptr));
}
// Number of occurrences of item in list
template<class T> int list_count (list<T>* lptr, const T& item)
{
  int count = 0;
  for (typename list<T>::const_iterator i=lptr->begin();
       i != lptr->end(); ++i) {
    if (*i == item) ++count;
  }
  return (count);
}
// Index of first occurrence of item in list
template<class T> int list_index (list<T>* lptr, const T& item)
{
  int ind = 0;
  typename list<T>::const_iterator i=lptr->begin();
  for (; i != lptr->end(); ++i, ++ind) {
    if (*i == item)
      break;
  }
  if (i == lptr->end())
    //_SWIG_exception(SWIG_ValueError, "index: item not in list");
    PyErr_SetString(PyExc_ValueError, "index: item not in list");
  return ind;
}
// Insert item at index.
// If index is negative insert at beginning
// If index is > end, insert at end
template<class T> void list_insert (list<T>* lptr, int index, const T& item)
{ 
  typename list<T>::iterator i;
  if (index <= 0)
    i = lptr->begin();
  else if (index >= lptr->size())
    i = lptr->end();
  else {
    int ind = 0;
    for (i = lptr->begin(); i != lptr->end(); ++i, ++ind) {
      if (ind == index)
        break;
    }
  }
  lptr->insert(i, item);
}
// Remove first occurrence of item
template<class T> void list_remove (list<T>* lptr, const T& item)
{
  int ind = 0;
  for (typename list<T>::iterator i=lptr->begin();
       i != lptr->end(); ++i, ++ind) {
    if (*i == item) {
      lptr->erase(i);
      return;
    }
  }
  //_SWIG_exception(SWIG_ValueError, "remove: item not in list");
    PyErr_SetString(PyExc_ValueError, "remove: item not in list");
}
// Get item at index
// Negative index looks from the end
template<class T> T list___getitem__ (list<T>* lptr, int index)
{
  T retval;
  int ind = 0;
  if (lptr->size() > 0) {
    if (index >= 0) {
      for (typename list<T>::iterator i=lptr->begin();
	   i != lptr->end(); ++i, ++ind) {
        if (ind == index) {
          retval = *i;
          break;
        }
      }
    }
    else {
      int ind = -1;
      for (typename list<T>::reverse_iterator i=lptr->rbegin();
	   i != lptr->rend(); ++i, --ind) {
        if (ind == index) {
          retval = *(i.base());
          break;
        }
      }
    }
  }
  if (ind != index)
    //_SWIG_exception(SWIG_IndexError, "__getitem__: index out of range");
    PyErr_SetString(PyExc_IndexError, "__getitem__: index out of range");
  return retval;
}
// Set item at index to value
// Negative index looks from the end
template<class T> void list___setitem__ (list<T>* lptr, int index, const T& value)
{
  if (lptr->size() > 0) {
    if (index >= 0) {
      int ind = 0;
      for (typename list<T>::iterator i=lptr->begin();
	   i != lptr->end(); ++i, ++ind) {
        if (ind == index) {
          typename list<T>::iterator j = lptr->insert(i, value);
          lptr->erase(++j);
          return;
        }
      }
    }
    else {
      int ind = -1;
      for (typename list<T>::reverse_iterator i=lptr->rbegin();
	   i != lptr->rend(); ++i, --ind) {
        if (ind == index) {
          typename list<T>::iterator j = lptr->insert(i.base(), value);
          lptr->erase(++j);
          return;
        }
      }
    }
  }
  //_SWIG_exception(SWIG_IndexError, "__setitem__: index out of range");
  PyErr_SetString(PyExc_IndexError, "__setitem__: index out of range");
}
// Delete item at index
// Negative index looks from the end
template<class T> void list___delitem__ (list<T>* lptr, int index)
{
  if (lptr->size() > 0) {
    if (index >= 0) {
      int ind = 0;
      for (typename list<T>::iterator i=lptr->begin();
	   i != lptr->end(); ++i, ++ind) {
        if (ind == index) {
          lptr->erase(i);
          return;
        }
      }
    }
    else {
      int ind = -1;
      for (typename list<T>::reverse_iterator i=lptr->rbegin();
	   i != lptr->rend(); ++i, --ind) {
        if (ind == index) {
          lptr->erase(i.base());
          return;
        }
      }
    }
  }
  //_SWIG_exception(SWIG_IndexError, "__delitem__: index out of range");
  PyErr_SetString(PyExc_IndexError, "__delitem__: index out of range");
}

/*
 * Functions to emulate a Python dict using a C++ map
 */

// Clear the dict
template<class T1, class T2> void dict_clear (map<T1, T2>* mptr)
{
  mptr->clear();
}
// Return a copy of the dict
template<class T1, class T2>  map<T1, T2> * dict_copy(map<T1, T2>* mptr)
{
  return new map<T1, T2> (mptr->begin(), mptr->end());
}
// Does the dict have this key?
template<class T1, class T2> int dict_has_key (map<T1, T2>* mptr, const T1& key)
{
  return (mptr->find(key) != mptr->end());
}
// Length of dict
template<class T1, class T2> int dict___len__ (map<T1, T2>* mptr)
{
  return mptr->size();
}
// Delete item with key
template<class T1, class T2> void dict___delitem__ (map<T1, T2>* mptr, const T1& key)
{
  typename map<T1, T2>::iterator i = mptr->find(key);
  if (i != mptr->end()) {
    mptr->erase(i);
  }
  else {
    //_SWIG_exception(SWIG_IndexError, "__delitem__: key not found");
    PyErr_SetString(PyExc_IndexError, "__delitem__: key not found");
  }
}
// Set an item
// This will replace existing items as Python does
template<class T1, class T2> void dict___setitem__ (map<T1, T2>* mptr, const T1& key, const T2& value)
{
  const typename map<T1, T2>::value_type v(key, value);
  pair<typename map<T1, T2>::iterator,bool> p = mptr->insert(v);
  if (!p.second) {
    mptr->erase(p.first);
    mptr->insert(v);
  }
}
// Get an item
template<class T1, class T2> T2 dict___getitem__ (map<T1, T2>* mptr, const T1& key)
{
  typename map<T1, T2>::iterator i = mptr->find(key);
  if (i != mptr->end()) {
    return (i->second);
  }
  //_SWIG_exception(SWIG_IndexError, "__getitem__: key not found");
  PyErr_SetString(PyExc_IndexError, "__getitem__: key not found");
}
// Adds (copies) all elements from dict d
// This will replace existing items as Python does
template<class T1, class T2> void dict_update (map<T1, T2>* mptr, const map<T1, T2>& d)
{
  for (typename map<T1, T2>::const_iterator ind = d.begin();
       ind != d.end(); ++ind) {
    const typename map<T1, T2>::value_type& v = *ind;
    pair<typename map<T1, T2>::iterator,bool> p = mptr->insert(v);
    if (!p.second) {
      mptr->erase(p.first);
      p = mptr->insert(v);
    }
  }
}
// This function will convert the map<T1, T2>
// to a list of (T1, T2) pairs, which later should be converted
// to a python list of tuples by a typemap
template<class T1, class T2> list<pair<T1, T2> > * dict_items(map<T1, T2>* mptr)
{
  list<pair<T1, T2> > *lp = new list<pair<T1, T2> >;
  for (typename map<T1, T2>::const_iterator ind = mptr->begin();
       ind != mptr->end(); ++ind)
  {
    lp->push_back((const pair<T1, T2>&) *ind);
  }
  return lp;
}
// This function will convert the map<T1, T2>
// to a list of (T1, T2*) pairs, which later need to be converted
// to a python list of tuples of (T1, T2) by a typemap
/*
template<class T1, class T2> list<pair<T1, T2*> > * dict_items(map<T1, T2>* mptr)
{
  list<pair<T1, T2*> > *lp = new list<pair<T1, T2*> >;
  for (map<T1, T2>::const_iterator ind = mptr->begin(); ind != mptr->end(); ++ind)
  {
    pair<T1, T2*> p(ind->first, (T2 *) &ind->second);
    lp->push_back(p);
  }
  return lp;
}
*/
// This function will return a list of T1 keys from map<T1, T2>,
// which later should be converted to a python list of keys by a typemap
template<class T1, class T2> list<T1> * dict_keys(map<T1, T2>* mptr)
{
  list<T1> *lk = new list<T1>;
  for (typename map<T1, T2>::const_iterator ind = mptr->begin();
       ind != mptr->end(); ++ind)
  {
    const T1& key = ind->first;
    lk->push_back(key);
  }
  return lk;
}
// This function will return a list of pointers to T2 values from map<T1, T2>,
// which later need to be converted to a python list of T2 references
// by an appropriate typemap
template<class T1, class T2> list<T2*> * dict_values(map<T1, T2>* mptr)
{
  list<T2*> *lv = new list<T2*>;
  for (typename map<T1, T2>::iterator ind = mptr->begin();
       ind != mptr->end(); ++ind)
  {
    T2& value = ind->second;
    lv->push_back(&value);
  }
  return lv;
}

#endif
