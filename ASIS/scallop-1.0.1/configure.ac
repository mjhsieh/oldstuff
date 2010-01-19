dnl -*- Autoconf -*-

dnl Process this file with autoconf to produce a configure script.
AC_INIT(src/exact_numerics.f)
AC_CONFIG_HEADERS(config.h)
AC_CONFIG_AUX_DIR(config-aux)
AC_CANONICAL_SYSTEM

dnl $Format: "AM_INIT_AUTOMAKE($Project$,$ProjectVersion$,[do-not-define])"$
AM_INIT_AUTOMAKE(scallop,1.0.1,[do-not-define])
dnl AM_CONFIG_HEADER(config.h)

dnl Checks for programs.
AC_PROG_RANLIB
AC_PROG_CXX(mpCC g++)
AC_PROG_F77(xlf f77 g77)
AC_F77_LIBRARY_LDFLAGS
AC_PROG_AWK

AM_CONDITIONAL(USING_GXX, test x$GXX = xyes)
AM_CONDITIONAL(USING_G77, test x$G77 = xyes)
AC_F77_WRAPPERS

if test x$GXX = xyes; then
   OTHER_CXXFLAGS="-frepo"
fi
AC_SUBST(OTHER_CXXFLAGS)

dnl AH_TEMPLATE([DEFINE_TEMPLATES],
dnl 	[This needs to be defined for g++ templates, at least the way KeLP does them.])
dnl if test $GXX = yes; then
dnl    AC_DEFINE(DEFINE_TEMPLATES)
dnl fi

dnl Checks for libraries.

AC_ARG_ENABLE(debug,
[  --enable-debug          Turn on debugging],
[case "${enableval}" in
  yes) debug=true ;;
  no)  debug=false ;;
  *) AC_MSG_ERROR(bad value ${enableval} for --enable-debug) ;;
esac],[debug=false])
AM_CONDITIONAL(DEBUG, test x$debug = xtrue)

if test x$debug = xtrue; then
   if test x$GXX = xyes; then
      CXXFLAGS="-g"
   else
      CXXFLAGS="-g"
   fi
   FFLAGS="-g"
else
   CXXFLAGS="-O2"
   FFLAGS="-O2"
fi

if test x$debug = xtrue; then
   if test x$G77 = xyes; then
      OTHER_FFLAGS="-fbounds-check -O2"
   fi
   if test x$F77 = xmpxlf || test x$F77 = xxlf; then
      OTHER_FFLAGS="-C"
   fi
else
   if test x$CXX = xmpCC || test x$CXX = xxlC; then
      OTHER_CXXFLAGS="-qarch=pwr3 -qtune=pwr3"
   fi
fi
AC_SUBST(OTHER_FFLAGS)

DEFAULT_FFTW_DIRS="/usr /usr/local"
DEFAULT_MPI_DIRS="/usr/local/mpi /usr/lpp/ppe.poe"
DEFAULT_MPICH_DIRS="/usr/lib/mpich /usr/local/mpich /scratch/slocal/gballs/extras/mpich"
DEFAULT_KELP_DIRS="$KELP_HOME /scratch/slocal/gballs/research/KeLP1.3_src /usr/local/apps/KeLP1.3"
dnl DEFAULT_HDF5_DIRS="/usr/local"

AC_ARG_WITH(fftw,
[  --with-fftw=PATH        Specify path to fftw installation],
fftw_dir="${withval}")
AC_ARG_WITH(kelp,
[  --with-kelp=PATH        Specify path to kelp installation],
kelp_dir="${withval}")
AC_ARG_WITH(mpi,
[  --with-mpi=PATH         Specify path to mpi installation],
mpi_dir="${withval}")
AC_ARG_WITH(mpich,
[  --with-mpich=PATH       Specify path to mpich installation],
mpich_dir="${withval}")
dnl AC_ARG_WITH(hdf5,
dnl [  --with-hdf5=PATH        Specify path to hdf5 installation],
dnl hdf5_dir="${withval}")

AC_MSG_CHECKING(for fftw library)
for DIRECTORY in ${fftw_dir} ${DEFAULT_FFTW_DIRS}; do
  if test -d $DIRECTORY; then
    if test -f $DIRECTORY/lib/librfftw.a && test -f $DIRECTORY/include/fftw.h; then
      AC_MSG_RESULT(found fftw library in $DIRECTORY)
      fftw_dir=$DIRECTORY
      FFTW_LD_DIR="-L${fftw_dir}/lib"
      FFTW_LD_LIBS="-lrfftw -lfftw"
      FFTW_INCLUDE_DIR="-I${fftw_dir}/include"
      AC_SUBST(FFTW_LD_DIR)
      AC_SUBST(FFTW_LD_LIBS)
      AC_SUBST(FFTW_INCLUDE_DIR)
      have_fftw_dir=true
      break
    fi
  fi
done
if test -z "$have_fftw_dir"; then
  AC_MSG_ERROR([cannot find the FFTW library...

Please re-run configure and set the path to your FFTW library installation,
e.g. /usr/local/lib, using the configure option --with-fftw=PATH])
fi

dnl I should actually check for the libraries, and that they work:
dnl at this point I only check that a known mpi or kelp directory exists.
AC_MSG_CHECKING(for mpi directory)
for DIRECTORY in ${mpi_dir} ${DEFAULT_MPI_DIRS}; do
  if test -d $DIRECTORY; then
    AC_MSG_RESULT(found $DIRECTORY)
    mpi_dir=$DIRECTORY
    MPI_LD_DIR="-L${mpi_dir}/lib"
    MPI_LD_LIBS="-lmpi"
    MPI_INCLUDE_DIR="-I${mpi_dir}/include"
    AC_SUBST(MPI_LD_DIR)
    AC_SUBST(MPI_LD_LIBS)
    AC_SUBST(MPI_INCLUDE_DIR)
    have_mpi_dir=true
    break
  fi
done
if test -z "$have_mpi_dir"; then
  AC_MSG_RESULT(no)
fi

if test -z "$have_mpi_dir"; then
  AC_MSG_CHECKING(for mpich directory)
  for DIRECTORY in ${mpich_dir} ${DEFAULT_MPICH_DIRS}; do
    if test -d $DIRECTORY; then
      AC_MSG_RESULT(found $DIRECTORY)
      mpich_dir=$DIRECTORY
      MPI_LD_DIR="-L${mpich_dir}/lib"
      MPI_LD_LIBS="-lmpich -lpmpich"
      MPI_INCLUDE_DIR="-I${mpich_dir}/include"
      AC_SUBST(MPI_LD_DIR)
      AC_SUBST(MPI_LD_LIBS)
      AC_SUBST(MPI_INCLUDE_DIR)
      have_mpich_dir=true
      break
    fi
  done
  if test -z "$have_mpich_dir"; then
    AC_MSG_RESULT(no)
  fi
fi

if test -z "${have_mpi_dir}" && test -z "${have_mpich_dir}"; then
  AC_MSG_ERROR([cannot find mpi library: 

Please re-run configure and set the path to your KeLP library installation,
e.g. /usr/lib/mpich or /usr/local/mpi, using the configure option
--with-mpi=PATH or --with-mpich=PATH])
fi

AC_MSG_CHECKING(for KeLP library directory)
for DIRECTORY in ${kelp_dir} ${DEFAULT_KELP_DIRS}; do
  if test -d $DIRECTORY; then
    AC_MSG_RESULT(found $DIRECTORY)
    kelp_dir=$DIRECTORY
    KELP_LD_DIR="-L${kelp_dir}/lib"
    KELP_LD_LIBS="-ldgrid -ldock -lkelp -larray -lmp++ -ltools -lz -lrt -lrfftw -lfftw"
    KELP_INCLUDE_DIR="-I${kelp_dir}/include"
    AC_SUBST(KELP_LD_DIR)
    AC_SUBST(KELP_LD_LIBS)
    AC_SUBST(KELP_INCLUDE_DIR)
    have_kelp_dir=true
    break
  fi
done

if test -z "${have_kelp_dir}"; then
  AC_MSG_ERROR([cannot find KeLP library...

Please re-run configure and set the path to your KeLP library installation,
e.g. /usr/local/KeLP1.3/lib, using the configure option --with-kelp=PATH])
fi

dnl AC_MSG_CHECKING(for HDF5 library)
dnl for DIRECTORY in $hdf5_dir, ${DEFAULT_HDF5_DIRS}; do
dnl   if test -d $DIRECTORY; then
dnl     if test -f $DIRECTORY/lib/libhdf5.a && test -f $DIRECTORY/include/hdf5.h
dnl       AC_MSG_RESULT(found hdf5 in $DIRECTORY)
dnl       hdf5_dir = $DIRECTORY
dnl       HDF5_LD_DIR="-L${hdf5_dir}/lib"
dnl       HDF5_LD_LIBS="-lhdf5"
dnl       HDF5_INCLUDE_DIR="-I${hdf5_dir}/include"
dnl       AC_SUBST(HDF5_LD_DIR)
dnl       AC_SUBST(HDF5_LD_LIBS)
dnl       AC_SUBST(HDF5_INCLUDE_DIR)
dnl       have_hdf5_dir=true
dnl       break
dnl     fi
dnl   fi
dnl done

dnl if test -z "${have_hdf5_dir}"; then
dnl   AC_MSG_ERROR([cannot find HDF5 library...

dnl Please re-run configure and set the path to your HDF5 library installation,
dnl e.g. /usr/local/hdf5, using the configure option --with-hdf5=PATH])
dnl fi

dnl Checks for header files.
AC_HEADER_STDC

dnl Checks for typedefs, structures, and compiler characteristics.
AC_C_CONST

dnl Checks for library functions.

dnl I should write a real check for this, but for now I'm just switching
dnl on machine type.
case "$target_os" in
  solaris*)
    NETWORKING_LIBS_DEFAULT="-lsocket -lnsl"
    ;;
  linux*)
    NETWORKING_LIBS_DEFAULT=""
    ;;
esac
NETWORKING_LIBS=${NETWORKING_LIBS_DEFAULT}
AC_SUBST(NETWORKING_LIBS)

AC_OUTPUT(Makefile kAPI/Makefile src/Makefile test/Makefile config-aux/Makefile)

