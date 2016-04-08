CONTENTS:
1: DESCRIPTION
  1.1: FEATURES
2: CONTENTS AND INSTALLATION
  2.1: OCTAVE
    2.1.1 - Install
    2.1.2 - Note for Windows users
  2.2: MATLAB
3: GETTING STARTED
  2.1: EXAMPLES
  2.2: HELP
4: HOW TO CONTRIBUTE
  4.1: CITING GeoPDEs
  4.2: PROVIDING FEEDBACK

-----------------------------------------------------

1. DESCRIPTION

  GeoPDEs is an open and free package for the research and teaching of Isogeometric Analysis, written in Octave and fully compatible with Matlab.

  The GeoPDEs package provides a common and flexible framework for implementing and testing new isogeometric methods for the solution of partial differential equations.

1.1 FEATURES

  These are some of the main features of GeoPDEs:

    * Dimension independent implementation: the same code is valid for curves, surfaces and volumes.
    * Implementation of multipatch domains, with conforming interfaces.
    * Div- and curl-conforming spline spaces, also in multipatch domains.
    * Examples for Poisson, linear elasticity, advection-diffusion, bilaplacian, Stokes and Maxwell equations.
    * Detailed documentation, accessible with the help command.

2. CONTENTS AND INSTALLATION

  Download and uncompress the file GeoPDEs_full.tar.gz. This file contains 
   all that is necessary to install GeoPDEs:
   - The GeoPDEs package: geopdes-3.0.1.tar.gz
   - The NURBS package: nurbs-1.3.12.tar.gz
   - The technical report: GeoPDES_report.pdf
   - This README file.
   - For Matlab users, the mex files of the NURBS toolbox: nurbs_mex_files.tar.gz

    
2.1. OCTAVE

2.1.1 - Install
 
 * Install the "nurbs" package by typing at the octave prompt
    pkg install nurbs-<version>.tar.gz

 * Install and load GeoPDEs by typing at the octave prompt
    pkg install geopdes-<version>.tar.gz
    pkg load geopdes

2.1.2 - Note for Windows users

   The mingw-based binary version of Octave distributed on
    http://octave.sf.net by default installs dynamically loaded binary
    functions in a different location as that of script functions. This
    can create problems when installing GeoPDEs. In this case, we
    recommend you to use the non-compiled version of the package. 

 * Install the nurbs package, like in step 1 of Section 1.1.1.

 * Uncompress and untar the GeoPDEs package, and add to the path the 
    folder "geopdes", including its subfolders, typing at the octave prompt:
    addpath (genpath ('geopdes'))

2.2. MATLAB

 * Uncompress and untar the NURBS and the GeoPDEs packages, in the files
    nurbs-<version>.tar.gz, geopdes-<version>.tar.gz, respectively.

 * Add the generated folders to the path, including their subfolders. 
    You can do this by typing in in the command window:
    addpath (genpath ('nurbs'))
    addpath (genpath ('geopdes'))

 * Install the mex-files for the "nurbs" package (OPTIONAL):
   - uncompress and untar the file nurbs_mex_files.tar.gz in the folder 
      "nurbs/inst" of the nurbs package
   - in Matlab, go to the folder "nurbs/inst" and run the script file 'compile'.
      This will compile the files and save the nurbs package to your Matlab path

3. GETTING STARTED

3.1. EXAMPLES

 * GeoPDEs contains a set of simple examples, that can be run with the commands:
    geopdes_base_examples
    geopdes_elasticity_examples
    geopdes_fluid_examples
    geopdes_maxwell_examples
 
    These will guide you through a set of menus to choose some examples for solving different problems. 

 * The source code and the data files for these examples can be found in the directories:
    geopdes/inst/examples
    geopdes/inst/solve

 * A collection of simple geometries constructed with NURBS can be found 
    in the following directory:

     geopdes/inst/examples/geometry_files

3.2. HELP

 * A detailed documentation for the software is not likely to be published.
    We suggest you to read the paper

      R. Vazquez.
      A new design for the implementation of isogeometric analysis
      in Octave and Matlab: GeoPDEs 3.0, Tech. report, IMATI-CNR, 2016

    In the paper we give an explanation of the architecture, the design and
    the main features of the code. This paper can be seen as brief user's guide.

 * All the functions in GeoPDEs contain a detailed description, that can be
    accessed with "help <function_name>". You can see a list of all the 
    functions in the package, typing in the Octave prompt:
     pkg describe -verbose geopdes
     help geopdes (if the package is not installed, also valid for Matlab users)

 * The format for the geometry files is explained in the files

     geopdes/doc/geo_specs_v21.txt      % Single patch geometry
     geopdes/doc/geo_specs_mp_v21.txt   % Multiple patches geometry

4. HOW TO CONTRIBUTE

4.1 CITING GeoPDEs

 * GeoPDEs has been developed as a part of our research. The best way to support the development of GeoPDEs is by citing one (or more) of the following references:

      R. Vazquez.
      A new design for the implementation of isogeometric analysis
      in Octave and Matlab: GeoPDEs 3.0, Tech. report, IMATI-CNR, 2016

      C. de Falco, A. Reali, R. Vazquez. 
      GeoPDES: a research tool for Isogeometric Analysis of PDEs, 
      Advances in Engineering Software, 42(12), 2011
      doi:10.1016/j.advengsoft.2011.06.010

    in any paper where GeoPDEs is used to obtain numerical results.

4.2 GIVING FEEDBACK

 * To give us feedback, report bugs, suggest new features or to ask to
    be directly involved in the further development of GeoPDEs, please
    subscribe to the mailing list following the instructions at
    https://lists.sourceforge.net/lists/listinfo/geopdes-users

 * Before posting a question to the mailing list, please browse the
    mailing list archives at
    https://sourceforge.net/mailarchive/forum.php?forum_name=geopdes-users
    to see whether it has already been answered before.
