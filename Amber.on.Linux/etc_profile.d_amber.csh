setenv AMBERHOME /opt/amber11
set path = ($path $AMBERHOME/exe)

source /opt/intel/Compiler/11.1/056/bin/iccvars.csh ia32
source /opt/intel/Compiler/11.1/056/bin/ifortvars.csh ia32
source /opt/intel/Compiler/11.1/056/mkl/tools/environment/mklvars32.csh
#source /opt/intel/cc/9.0/bin/iccvars.csh
#source /opt/intel/fc/9.0/bin/ifortvars.csh
#setenv MKL_HOME /opt/intel/mkl721
setenv LD_LIBRARY_PATH ${LD_LIBRARY_PATH}:${MKL_HOME}/lib/32
