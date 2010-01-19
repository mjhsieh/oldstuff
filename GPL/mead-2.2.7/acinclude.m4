dnl Autoconf macros used for MEAD configuration.
dnl ID: $Id: acinclude.m4,v 1.14 2005/11/10 16:05:41 bashford Exp $
dnl SOURCE: $Source: /cvs-repository/bashford/cvsroot/mead/acinclude.m4,v $

dnl --------------------------------------------------------------------
dnl Python


AC_DEFUN([AC_PYTHON_NUMPY],
  [ AC_MSG_CHECKING(for python and Numeric API)
    AC_MSG_RESULT([])
    ac_python_numpy_p=no
    for pytry in python2.4 python2.3 python2.1 python2 python; do
      AC_CHECK_PROGS(PYTHON_BIN, [$pytry])
      if test x$PYTHON_BIN != x; then
        ac_python_prefix=`$pytry -c 'import sys; print sys.prefix'`
        [ac_python_version=`$pytry -c 'import sys; print sys.version[:3]'`]
        ac_python_incdir=$ac_python_prefix/include/python$ac_python_version
        AC_CHECK_HEADER([$ac_python_incdir/Python.h], 
                        [ac_python_h_found=yes],
                        [ac_python_h_found=no])
       
        if test x$ac_python_h_found = xyes; then
          AC_CHECK_HEADER([$ac_python_incdir/Numeric/arrayobject.h], 
                          [ac_numeric_arrayobject_h_found=yes],
                          [ac_numeric_arrayobject_h_found=no],
                          [[
                            #include <$ac_python_incdir/Python.h>
                          ]] )
          if test x$ac_numeric_arrayobject_h_found = xyes; then
            ac_python_numpy_p=yes
            AC_ARG_WITH(py-site-packages-dir,
              [  --with-py-site-packages-dir   installation dir for python packages],
              [PYSITEPKGDIR="$withval"],
              [PYSITEPKGDIR=])

            if test x$PYSITEPKGDIR = x ; then
              ac_python_exec_prefix=`$pytry -c 'import sys; print sys.exec_prefix'`
              PYSITEPKGDIR=$ac_python_prefix/lib/python$ac_python_version/site-packages
              if test -d $PYSITEPKGDIR; then
                true
              else
                PYSITEPKGDIR=""
                AC_MSG_WARN([Failed to find python site-packages dir])
                AC_MSG_WARN([If you want to install PyMead pyton package,])
  	      AC_MSG_WARN([specify location using --with-py-site-packages-dir])
              fi
            fi
            break
          fi
        fi
      fi
    done
    AC_MSG_RESULT([Results of Python and Numeric API check:])
    if test x$ac_python_numpy_p = xyes; then
      PYTHONINC=$ac_python_incdir
      AC_MSG_RESULT([  Python binary having Numeric:  $PYTHON_BIN])
      AC_MSG_RESULT([  Python include dir:            $PYTHONINC])
      AC_SUBST(PYTHONINC)
      AC_SUBST(PYSITEPKGDIR)
    else
      AC_MSG_WARN([   No Python with Numeric and API found])
      AC_MSG_WARN([   Python interface will not be built])
    fi
  ])  
      




dnl --------------------------------------------------------------------
dnl Guile

AC_DEFUN([AC_PROG_GUILE], [
  AC_PATH_PROG(GUILE, guile)
  if test -z "$GUILE" ; then
    AC_MSG_ERROR([Guile interpreter not found, \$GUILE not set])
  fi
])

AC_DEFUN([AC_GUILE_DIRS], [
  AC_PROG_GUILE
])

dnl --------------------------------------------------------------------
dnl Wrapping

AC_DEFUN([AC_PROG_SWIG], [
  AC_PATH_PROG(SWIG, swig)
  if test -z "$SWIG" ; then
    echo WARNING: SWIG not found, \$SWIG not set
    SWIG=
  fi
])

AC_DEFUN([AC_PROG_PERL], [
  AC_PATH_PROG(PERL, perl)
  if test -z "$PERL" ; then
    echo WARNING: PERL not found, \$PERL not set
    PERL=
  fi
])

AC_DEFUN([AC_PROG_PATCH], [
  AC_PATH_PROG(PATCH, patch)
  if test -z "$PATCH" ; then
    echo WARNING: PATCH not found, \$PATCH not set
    PATCH=
  fi
])

dnl acinclude.m4 ends here
