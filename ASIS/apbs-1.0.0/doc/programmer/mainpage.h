/** 
 * <center>
 * @mainpage APBS Programmers Guide
 *  
 * APBS was written by Nathan A. Baker.<br>
 * Additional contributing authors listed in the code documentation.
 * </center>
 * 
 * <hr width="100%">
 * @section toc   Table of Contents
 * <ul> 
 *   <li> @ref style
 *   <li> @ref api
 *     <ul>
 *       <li> <a href="modules.html">Modules</a>
 *       <li> <a href="annotated.html">Class list</a>
 *       <li> <a href="functions.html">Class members</a>
 *       <li> <a href="globals.html">Class methods</a>
 *    </ul>
 * </ul>
 * 
 * <hr width="100%">
 * @section license License
 *
 * Primary author: <a href="http://agave.wustl.edu/">Nathan A. Baker</a> (<a href="mailto:baker@biochem.wustl.edu">baker@biochem.wustl.edu</a>)<br>
 * Department of Biochemistry and Molecular Biophysics<br>
 * Center for Computational Biology<br>
 * Washington University in St. Louis<br>
 * Additional contributing authors are listed in the code documentation.<br>
 * 
 * <p> Copyright (c) 2002-2008, Washington University in St. Louis.
 * Portions Copyright (c) 2002-2008.  Nathan A. Baker
 * Portions Copyright (c) 1999-2002.  The Regents of the University of California.
 * Portions Copyright (c) 1995.  Michael Holst
 * 
 * <p> All rights reserved.
 * 
 * <p> Redistribution and use in source and binary forms, with or without
 * modification, are permitted provided that the following conditions are met: 
 * 
 * <ul>
 * <li> Redistributions of source code must retain the above copyright notice, this
 * list of conditions and the following disclaimer.  
 * 
 * <li> Redistributions in binary form must reproduce the above copyright notice,
 * this list of conditions and the following disclaimer in the documentation
 * and/or other materials provided with the distribution.
 * 
 * <li> Neither the name of Washington University in St. Louis nor the names of its
 * contributors may be used to endorse or promote products derived from this
 * software without specific prior written permission.
 * </ul>
 * <p>
 * THIS SOFTWARE IS PROVIDED BY THE COPYRIGHT HOLDERS AND CONTRIBUTORS
 * "AS IS" AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING, BUT NOT
 * LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY AND FITNESS FOR
 * A PARTICULAR PURPOSE ARE DISCLAIMED. IN NO EVENT SHALL THE COPYRIGHT OWNER OR
 * CONTRIBUTORS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL, SPECIAL,
 * EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED TO,
 * PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
 * PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
 * LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
 * NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
 * SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
 * 
 * <hr>
 * <p>
 * This documentation provides information about the programming interface
 * provided by the APBS software and a general guide to linking to the APBS
 * libraries.  Information about installation, configuration, and general usage
 * can be found in the <a href="user.html">User's Guide</a>.
 * <hr>
 * 
 *  @section style Programming Style
 * 
 *  <p>
 *  APBS was developed following the <a
 *  href="http://scicomp.ucsd.edu/~mholst/codes/maloc/cleanc.html">Clean OO
 *  C</a> style of Mike Holst.  In short, Clean OO C code is written in a
 *  object-oriented, ISO C-compliant fashion, and can be compiled with either a
 *  C or C++ compiler.  <p> Following this formalism, all public data is
 *  enclosed in structures which resemble C++ classes.  These structures and
 *  member functions are then declared in a public header file which provides a
 *  concise description of the interface for the class.  Private functions and
 *  data are included in private header files (or simply the source code files
 *  themselves) which are not distributed.  When using the library, the
 *  end-user only sees the public header file and the compiled library and is
 *  therefore (hopefully) oblivious to the private members and functions.  Each
 *  class is also equipped with a constructor and destructor function which is
 *  responsible for allocating and freeing any memory required by the
 *  instatiated objects.
 * 
 *  <p>
 *  As mentioned above, public data members are enclosed in C structures which
 *  are visible to the end-user.  Public member functions are generated by
 *  mangling the class and function names <i>and</i> passing a pointer to the
 *  object on which the member function is supposed to act.  For example, a
 *  public member function with the C++ declaration 
 *    <pre>
 *   public double Foo::bar(int i, double d)
 *   </pre>
 * would be declared as
 *   <pre>
 *   VEXTERNC double Foo_bar(Foo *thee, int i, double d)
 *   </pre>
 * where <code>VEXTERNC</code> is a compiler-dependent macro, the underscore
 * <code>_</code> replaces the C++ double-colon <code>::</code>, and
 * <code>thee</code> replaces the <code>this</code> variable implicit in all
 * C++ classes.  Since they do not appear in public header files, private
 * functions could be declared in any format pleasing to the user, however, the
 * above declaration convention should generally be used for both public and
 * private functions.  Within the source code, the public and private function
 * declarations/definitions are prefaced by the macros <code>VPUBLIC</code> and
 * <code>VPRIVATE</code>, respectively.  These are macros which reduce global
 * name pollution, similar to encapsulating private data withing C++ classes.
 *  
 * <p>
 * The only C++ functions not explicitly covered by the above declaration
 * scheme are the constructors (used to allocate and initialize class data
 * members) and destructors (used to free allocated memory).  These are
 * declared in the following fashion:  a constructor with the C++ declaration
 *    <pre>
 *    public void Foo::Foo(int i, double d)
 *    </pre>
 * would be declared as
 *     <pre>
 *     VEXTERNC Foo* Foo_ctor(int i, double d)
 *     </pre>
 * which returns a pointer to the newly constructed <code>Foo</code> object.
 * Likewise, a destructor declared as
 *     <pre>
 *     public void Foo::~Foo()
 *     </pre>
 * in C++ would be
 *     <pre>
 *     VEXTERNC void Foo_dtor(Foo **thee)
 *     </pre>
 * in Clean OO C.
 * <p>
 * Finally, inline functions in C++ are simply treated as macros in Clean OO C
 * and declared/defined using <code>#define</code> statements in the public
 * header file.
 * <p>
 * See any of the APBS header files for more information on Clean OO C
 * programming styles.
 * 
 * @section api Application programming interface documentation
 * <p>
 * The API documentation for this code was generated by <a
 * href="http://www.doxygen.org">doxygen</a>.  You can either view the API
 * documentation by using the links at the top of this page, or the slight
 * re-worded/re-interpreted list below:
 *    <ul>
 *    <li> <a href="modules.html">Class overview</a>
 *    <li> <a href="annotated.html">Class declarations</a>
 *    <li> <a href="functions.html">Class members</a>
 *    <li> <a href="globals.html">Class methods</a>
 *    </ul>
 * 
 */
