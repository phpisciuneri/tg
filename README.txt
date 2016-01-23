Hints:
=======
To avoid "No errors found" notifications, declare BOOST_PRG_MON_CONFIRM=no environment variable.


Build Instructions for xenon and xenon01
========================================
    
0. Load the proper compilers etc.  I have put them all in modules and you can load them all at once with:
   * "module load prog/intel" which essentially loads cmake, boost, intel compilers, mpi, lapack etc.
	 
1. Check out mytools (we will link with this package) in some directory where you want your source code to be:
   * "svn co https://collab.sam.pitt.edu/svn/mytools/trunk mytools"
	    
2. Build mytools:
   * create a build directory for mytools
   * from your mytools build directory run: "cmake -C ~/path/to/mytools/site-macros/caches/intel.cmake ~/path/to/mytools/"
   * this should be successful, so type "make" to build
	     
3. Check out the IPLMCFD source code in some directory where you want your source code to be:
   * "svn co https://collab.sam.pitt.edu/svn/iplmcfd/trunk iplmcfd"
		
4. Build iplmcfd:
   * create a build directory for iplmcfd
   * from within your iplmcfd build directory run: "cmake -C ~/path/to/iplmcfd/site-macros/caches/intel.cmake ~/path/to/iplmcfd/"
   * cmake will complain about not finding mytools
   * now run "ccmake ."
   * change the "mytools_DIR" to the location of your mytools build from step 2
   * type "c" to configure
   * turn "IPLMCFD_IMPLICIT_MPI" to "OFF"
   * turn "ODEPACK_IMPLICIT_BLAS" to "OFF"
   * type "c" to configure
   * type "c" to configure again
   * now type "g" to generate
   * the interface will close and you can run "make"
				    
5. Build should be successful
							    
6. ... More to come on how to run the code soon ... 



Build Instructions for Visual Studio

====================================
0. Tested with: CMake version: 2.8.5, Visual Studio 2010, Intel Fortran XE 12.0, MPICH2 1.4.1 Windows build

1. Build mytools
	* Open cmake-gui
	* Point to source and build directories
	* Configure - Select "Visual Studio 10" and select "Use default native compilers"	
	* Set the following options under the CMAKE group ( see "\Modules\Platform\Windows-Intel-Fortran.cmake" of CMake install directory ):
     (Last argument is mixed C++/Fortran call requirement peculiar to iplmcfd in windoze - until fixed)
         CMAKE_Fortran_FLAGS                /W1 /nologo /fpp /libs:dll /threads /iface:mixed_str_len_arg
         CMAKE_Fortran_FLAGS_DEBUG          /debug:full /dbglibs
         CMAKE_Fortran_FLAGS_MINSIZEREL     /O2 /D NDEBUG
         CMAKE_Fortran_FLAGS_RELEASE        /O1 /D NDEBUG
         CMAKE_Fortran_FLAGS_RELWITHDEBINFO /O1 /debug:full /D NDEBUG
    * Configure & generate
	* Open "mytools.sln" in build directory to build solution in VS10
	
2. Build IPLMCFD
	* Open cmake-gui
	* Point to source and build directions
	* Configure - Select "Visual Studio 10" and select "Use default native compilers"
	* Set the following options under the CMAKE group ( see "\Modules\Platform\Windows-Intel-Fortran.cmake" of CMake install directory ):
		CMAKE_Fortran_FLAGS                /W1 /nologo /fpp /libs:dll /threads /iface:mixed_str_len_arg
		CMAKE_Fortran_FLAGS_DEBUG          /debug:full /dbglibs
		CMAKE_Fortran_FLAGS_MINSIZEREL     /O2 /D NDEBUG
		CMAKE_Fortran_FLAGS_RELEASE        /O1 /D NDEBUG
		CMAKE_Fortran_FLAGS_RELWITHDEBINFO /O1 /debug:full /D NDEBUG
	*  Set TVMET_INCLUDE_DIRS under Ungrouped Entries
	*  Set Subversion_SVN_EXECUTABLE under Ungrouped Entries
	*  Deselect IPLMCFD_IMPLICIT_MPI under IPLMCFD group
	*  Remove cxx.lib from MPI_CXX_LIBRARIES under MPI group
	*  Deselect ODEPACK_IMPLICIT_BLAS under Ungrouped Entries
	*  Set BLAS_mkl_LIBRARY under BLAS group to point to mkl_rt.lib of intel installation - use ia32
	*  Set BLAS_guide_LIBRARY under Ungrouped Entries to point to mkl_rt.lib of intel installation - use ia32
	*  Set LAPACK_mkl_lapack_LIBRARY under LAPACK group to point to mkl_rt.lib of intel installation - use ia32
 	*  Configure & Generate
    *  Open "iplmcfd.dirty.sln" in build directory to build solution in VS10
	*  If there is the following link error: "LINK : fatal error LNK1104: cannot open file 'libboost_filesystem-vc100-mt-gd-1_46_1.lib'"
		Go to Property Manager - iplmcfd - Debug | Win32 - Microsoft.Cpp.Win32.user
		Right click and select Properties
		Under Additional Library Directories add Boost lib directory
	*  If there is the following link error: "LINK : fatal error LNK1104: cannot open file 'ifconsol.lib'"
		Go to Property Manager - iplmcfd - Debug | Win32 - Microsoft.Cpp.Win32.user
		Right click and select Properties
		Under Additional Library Directories add intel compiler lib directory



CommitTicketUpdater -- Update tickets based on commit messages.
==============================================================
(copied from Trac plugin docs)

This component hooks into changeset notifications and searches commit messages for text in the form of:

command #1
command #1, #2
command #1 & #2 
command #1 and #2
Instead of the short-hand syntax "#1", "ticket:1" can be used as well, e.g.:

command ticket:1
command ticket:1, ticket:2
command ticket:1 & ticket:2 
command ticket:1 and ticket:2
In addition, the ':' character can be omitted and issue or bug can be used instead of ticket.

You can have more than one command in a message. The following commands are supported. There is more than one spelling for each command, to make this as user-friendly as possible.

close, closed, closes, fix, fixed, fixes
The specified tickets are closed, and the commit message is added to them as a comment.
references, refs, addresses, re, see
The specified tickets are left in their current status, and the commit message is added to them as a comment.
A fairly complicated example of what you can do is with a commit message of:

eg: Changed blah and foo to do this or that. Fixes #10 and #12, and refs #12.

This will close #10 and #12, and add a note to #12.


Using SVNMERGE tool for merges
==============================

SVNMERGE is a neat tool to automate merging, making it more reasonable. 
See quick start tutorial at http://www.orcaware.com/svn/wiki/Svnmerge.py. 

SVNMERGE is available on linuxed (just some python script, download from the web or
just use distro installer). Binaries for Windoze are also available on the website, 
but they don't seem to work that cleanly. Eitherway, it should be ok to do any kind of 
merging activity on linux. 

In a nutshell

1. Normal workflow: create a feature branch
2. Checkout a fresh copy of the feature branch and run "svnmerge init" on it.
3. The tool just modifies properties of the top level dir, it doesn't create or edit any files or anything. 
   Commit after init with new property. 
4. Pull to trunk or from trunk to branches and so on based on the normal workflow. See tutorial online./ 


