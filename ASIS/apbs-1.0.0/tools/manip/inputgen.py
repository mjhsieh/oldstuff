""" inputgen class

    Create an APBS input file using psize data

    Written by Todd Dolinsky based on original sed script by Nathan Baker

        ----------------------------

    Version:  $Id: inputgen.py 1251 2008-03-30 22:32:00Z sobolevnrm $
   
	APBS -- Adaptive Poisson-Boltzmann Solver

	  Nathan A. Baker (baker@biochem.wustl.edu)
	  Dept. Biochemistry and Molecular Biophysics
	  Center for Computational Biology
	  Washington University in St. Louis

	  Additional contributing authors listed in the code documentation.

	Copyright (c) 2002-2008, Washington University in St. Louis.
	Portions Copyright (c) 2002-2008.  Nathan A. Baker
	Portions Copyright (c) 1999-2002.  The Regents of the University of California.
	Portions Copyright (c) 1995.  Michael Holst

	All rights reserved.

	Redistribution and use in source and binary forms, with or without
	modification, are permitted provided that the following conditions are met: 

	* Redistributions of source code must retain the above copyright notice, this
	list of conditions and the following disclaimer.  

	* Redistributions in binary form must reproduce the above copyright notice,
	this list of conditions and the following disclaimer in the documentation
	and/or other materials provided with the distribution.

	* Neither the name of Washington University in St. Louis nor the names of its
	contributors may be used to endorse or promote products derived from this
	software without specific prior written permission.

	THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
	"AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
	LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
	A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
	CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
	EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
	PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
	PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
	LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
	NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
	SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.


    ----------------------------
"""

# User - Definable Variables: Default values

# cfac = 1.7                  # Factor by which to expand mol dims to
                              # get coarse grid dims
# fadd = 20                   # Amount to add to mol dims to get fine
                              # grid dims
# space = 0.50                # Desired fine mesh resolution
# gmemfac = 200               # Number of bytes per grid point required 
                              # for sequential MG calculation 
# gmemceil = 400              # Max MB allowed for sequential MG
                              # calculation.  Adjust this to force the
                              # script to perform faster calculations (which
                              # require more parallelism).
# ofrac = 0.1                  # Overlap factor between mesh partitions
# redfac = 0.25               # The maximum factor by which a domain
                              # dimension can be reduced during focusing

import string, sys
import psize


class inputGen:
    """
        DEPRECATED:  Please use the Input class instead!  See the
                     main() function for usage.
    """
    def __init__(self, pqrpath, size, method, async):
        size.runPsize(pqrpath)
        self.input = Input(pqrpath, size, method, async)

    def printInput(self):
        self.input.printInputFiles()

class Elec:
    """
        An object for the ELEC section of an APBS input file
    """
    def __init__(self, size, method, asyncflag):
        """
            Initialize the variables that can be set in this object
            Users can modify any of these variables (that's why
            they're here!)
        """

        # If this is an async or parallel calc, we want to use
        # the per-grid dime rather than the global dime.
        
        self.dime = size.getFineGridPoints()
        gmem = 200.0 * self.dime[0] * self.dime[1] * self.dime[2] / 1024.0 / 1024.0
        if method == "": # method not named - use ceiling
            if gmem > size.getConstant("gmemceil"): method = "mg-para"
            else: method = "mg-auto"

        if method == "mg-para":
            self.dime = size.getSmallest()

        self.method = method
        self.glen = size.getCoarseGridDims()
        self.cglen = size.getCoarseGridDims()
        self.fglen = size.getFineGridDims()
        self.pdime = size.getProcGrid()
        
        self.label = ""
        self.nlev = 4
        self.ofrac = 0.1
        self.async = 0
        self.asyncflag = asyncflag
        self.cgcent = "mol 1"
        self.fgcent = "mol 1"
        self.gcent = "mol 1"
        self.mol = 1
        self.lpbe = 1
        self.npbe = 0
        self.bcfl = "sdh"
        self.ion = [] # Multiple ions possible
        self.pdie = 2.0
        self.sdie = 78.54
        self.srfm = "smol"
        self.chgm = "spl2"
        self.sdens = 10.0
        self.srad = 1.4
        self.swin = 0.3
        self.temp = 298.15
        self.gamma = 0.105
        self.calcenergy = "total"
        self.calcforce = "no"
        self.write = [] # Multiple write statements possible
    
    def __str__(self):
        """
            Return the elec statement as a string. Check the method
            to see which keywords to use.
        """
        text = "elec %s\n" % self.label
        text += "    %s\n" % self.method
        text += "    dime %i %i %i\n" % (self.dime[0], self.dime[1], self.dime[2])
        if self.method == "mg-manual":
            text += "    nlev %i\n" % self.nlev
            text += "    glen %.3f %.3f %.3f\n" % (self.glen[0], self.glen[1], self.glen[2])
            text += "    gcent %s\n" % self.gcent
        elif self.method == "mg-auto":
            text += "    cglen %.4f %.4f %.4f\n" % (self.cglen[0], self.cglen[1], self.cglen[2])
            text += "    fglen %.4f %.4f %.4f\n" % (self.fglen[0], self.fglen[1], self.fglen[2])
            text += "    cgcent %s\n" % self.cgcent
            text += "    fgcent %s\n" % self.fgcent
        elif self.method == "mg-para":
            text += "    pdime %i %i %i\n" % (self.pdime[0], self.pdime[1], self.pdime[2])
            text += "    ofrac %.1f\n" % self.ofrac
            text += "    cglen %.4f %.4f %.4f\n" % (self.cglen[0], self.cglen[1], self.cglen[2])
            text += "    fglen %.4f %.4f %.4f\n" % (self.fglen[0], self.fglen[1], self.fglen[2])
            text += "    cgcent %s\n" % self.cgcent
            text += "    fgcent %s\n" % self.fgcent
            if self.asyncflag == 1:
                text += "    async %i\n" % self.async
        text += "    mol %i\n" % self.mol
        if self.lpbe: text += "    lpbe\n"
        else: text += "    npbe\n"
        text += "    bcfl %s\n" % self.bcfl
        for ion in self.ion:
            text += "    ion %.2f %.3f %.2f\n" % (ion[0], ion[1], ion[2])               
        text += "    pdie %.4f\n" % self.pdie                
        text += "    sdie %.4f\n" % self.sdie                
        text += "    srfm %s\n" % self.srfm                   
        text += "    chgm %s\n" % self.chgm
        text += "    sdens %.2f\n" % self.sdens
        text += "    srad %.2f\n" % self.srad          
        text += "    swin %.2f\n" % self.swin         
        text += "    temp %.2f\n" % self.temp     
        text += "    gamma %.3f\n" % self.gamma    
        text += "    calcenergy %s\n" % self.calcenergy
        text += "    calcforce %s\n" % self.calcforce
        for write in self.write:
            text += "    write %s %s %s\n" % (write[0], write[1], write[2])
        text += "end\n"
        return text
        
class Input:
    """
        The input class.  Each input object is one APBS input file.
    """

    def __init__(self, pqrpath, size, method, asyncflag):
        """
            Initialize the input file class.  Each input file contains
            a PQR name, a list of elec objects, and a list of strings
            containing print statements.  For starters assume two
            ELEC statements are needed, one for the inhomgenous and
            the other for the homogenous dielectric calculations.

            Users can edit the elec statements and the print statements.

            This assumes you have already run psize, either by
                 size.runPsize(/path/to/pqr) or

                 size.parseString(string)
                 size.setAll()

            Parameters
                pqrpath:   The path to the PQR file (string)
                size:      The Psize object (psize)
                method:    The method (para, auto, manual, async) to use
                asyncflag: 1 if async is desired, 0 otherwise
        """ 

        self.pqrpath = pqrpath
        self.asyncflag = asyncflag

        # Initialize variables to default elec values

        elec1 = Elec(size, method, asyncflag)
        elec2 = Elec(size, method, asyncflag)
        setattr(elec2, "sdie", 2.0)
        self.elecs = [elec1, elec2]
     
        i = string.rfind(pqrpath, "/") + 1
        self.pqrname = pqrpath[i:]

        self.prints = ["print energy 2 - 1 end"]     

    def __str__(self):
        """
            Return the text of the input file
        """
        text  = "read\n"
        text += "    mol pqr %s\n" % self.pqrname
        text += "end\n"
        for elec in self.elecs:
            text += str(elec)            
        for prints in self.prints:
            text += prints
        text += "\nquit\n"
        return text
  
    def printInputFiles(self):
        """
            Make the input file(s) associated with this object
        """
        period = string.find(self.pqrpath,".")
        if self.asyncflag == 1:
            outname = self.pqrpath[0:period] + "-para.in"

            # Temporarily disable async flag
            for elec in self.elecs:
                elec.asyncflag = 0
            file = open(outname, "w")
            file.write(str(self))
            file.close()

            # Now make the async files
            elec = self.elecs[0]
            
            nproc = elec.pdime[0] * elec.pdime[1] * elec.pdime[2]
            for i in range(int(nproc)):
                outname = self.pqrpath[0:period] + "-PE%i.in" % i
                for elec in self.elecs:
                    elec.asyncflag = 1
                    elec.async = i
                file = open(outname, "w")
                file.write(str(self))
                file.close()
        
        else:
            if period > 0:
                outname = self.pqrpath[0:period] + ".in"
            else:
                outname = self.pqrpath + ".in"
            file = open(outname, "w")
            file.write(str(self))
            file.close()

def splitInput(filename):
    """
        Split the parallel input file into multiple async file names

        Parameters
            filename:  The path to the original parallel input
                       file (string)
    """
    nproc = 0
    file = open(filename)
    text = ""
    while 1:
        line = file.readline()
        if line == "": break
        text += line
        line = string.strip(line)
        if line.startswith("pdime"): # Get # Procs
            words = string.split(line)
            nproc = int(words[1]) * int(words[2]) * int(words[3])

    if nproc == 0:
        sys.stderr.write("%s is not a valid APBS parallel input file!\n" % filename)
        sys.stderr.write("The inputgen script was unable to asynchronize this file!\n")
        sys.exit(2)

    period = string.find(filename,".")
    for i in range(nproc):
        outname = filename[0:period] + "-PE%i.in" % i
        outtext = string.replace(text, "mg-para\n","mg-para\n    async %i\n" % i)
        outfile = open(outname, "w")
        outfile.write(outtext)
        outfile.close()
          
def usage():
    """
        Display the usage information for this script
    """
    size = psize.Psize()
    usage = "\n"
    usage = usage + "Use this script to generate new APBS input files or split an existing\n"
    usage = usage + "parallel input file into multiple async files.\n\n"
    usage = usage + "Usage: inputgen.py [opts] <filename>\n"
    usage = usage + "Optional Arguments:\n"
    usage = usage + "  --help               : Display this text\n"
    usage = usage + "  --split              : Split an existing parallel input file to multiple\n"
    usage = usage + "                         async input files.\n"
    usage = usage + "  --method=<value>     : Force output file to write a specific APBS ELEC\n"
    usage = usage + "                         method.  Options are para (parallel), auto\n"
    usage = usage + "                         (automatic), manual (manual), or async (asynchronous).\n"
    usage = usage + "                         solve.  async will result in multiple input files.\n"
    usage = usage + "  --cfac=<value>       : Factor by which to expand molecular dimensions to\n"
    usage = usage + "                         get coarse grid dimensions.\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("cfac")
    usage = usage + "  --fadd=<value>       : Amount to add to molecular dimensions to get fine\n"
    usage = usage + "                         grid dimensions.\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("fadd")
    usage = usage + "  --space=<value>      : Desired fine mesh resolution\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("space")
    usage = usage + "  --gmemfac=<value>    : Number of bytes per grid point required\n"
    usage = usage + "                         for sequential MG calculation\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("gmemfac")
    usage = usage + "  --gmemceil=<value>   : Max MB allowed for sequential MG\n"
    usage = usage + "                         calculation.  Adjust this to force the\n"
    usage = usage + "                         script to perform faster calculations (which\n"
    usage = usage + "                         require more parallelism).\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("gmemceil")
    usage = usage + "  --ofrac=<value>      : Overlap factor between mesh partitions\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("ofrac")
    usage = usage + "  --redfac=<value>     : The maximum factor by which a domain\n"
    usage = usage + "                         dimension can be reduced during focusing\n"
    usage = usage + "                         [default = %g]\n" % size.getConstant("redfac")
    sys.stderr.write(usage)
    sys.exit(2)

def main():

    import getopt
    filename = ""
    shortOptList = ""
    longOptList = ["help","split","method=","cfac=","space=","gmemceil=","gmemfac=","ofrac=","redfac="]

    try:
        opts, args = getopt.getopt(sys.argv[1:], shortOptList, longOptList)
    except getopt.GetoptError, details:
        sys.stderr.write("Option error (%s)!\n" % details)
        usage()
        
    if len(args) != 1:
        sys.stderr.write("Invalid argument list!\n")
        usage()
    else:
        filename = args[0]

    method = ""
    size = psize.Psize()
    async = 0
    split = 0
    
    for o, a in opts:
        if o == "--help":
            usage()
        if o == "--split": split = 1
        if o == "--method":
            if a == "para":
                sys.stdout.write("Forcing a parallel calculation\n")
                method = "mg-para"
            elif a == "auto":
                sys.stdout.write("Forcing a sequential calculation\n")
                method = "mg-auto"
            elif a == "async":
                sys.stdout.write("Forcing an asynchronous calculation\n")
                method = "mg-para"
                async = 1
            elif a == "manual":
                sys.stdout.write("Forcing a manual calculation\n")
                method = "mg-manual"
            else:
                sys.stdout.write("Incorrect method argument: %s\n" % a)
                sys.stdout.write("Defaulting to memory dependent result\n")
        if o == "--cfac":
            size.setConstant("cfac", float(a))
        if o == "--space":
            size.setConstant("space", float(a))
        if o == "--gmemfac":
            size.setConstant("gmemfac", int(a))
        if o == "--gmemceil":
            size.setConstant("gmemceil",  int(a))
        if o == "--ofrac":
            size.setConstant("ofrac", float(a))
        if o == "--redfac":
            size.setConstant("redfac", float(a))

    if split == 1:
        splitInput(filename)
    else:
        size.runPsize(filename)
        input = Input(filename, size, method, async)
        input.printInputFiles()

if __name__ == "__main__": main()