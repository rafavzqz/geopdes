CONTENTS:
0: CHANGELOG
1: INSTALLATION
  1.1: OCTAVE
    1.1.1 - Install
    1.1.2 - Test
    1.1.3 - Note for windows users
  1.2: MATLAB
2: GETTING STARTED
  2.1: EXAMPLES
  2.2: HELP
3: HOW TO CONTRIBUTE
  3.1: CITING GeoPDEs
  3.2: PROVIDING FEEDBACK

-----------------------------------------------------

0: CHANGELOG

Below is a list of the main changes introduced with each new release

Version geopdes_***-1.1.0

* Changed the way the examples are called.
* Make use of tensor product in msh_push_forward_2d(3d) and geo_2d(3d)_nurbs.
    This requires version 1.3.4 of the nurbs toolbox.
* New file format for multipatch geometries (v.0.7).
* Changed the way the multipatch numbering (and orientation) is set.
* Added a new function to compute the deformed geometry in geopdes_elasticity.
* Added multipatch examples in geopdes_elasticity.
* Added Subgrid method in geopdes_fluid.
* Added multipatch examples for TH and SG elements in geopdes_fluid.
* Modified the C version of the operators to be compatible with Octave 3.4.
* Modified the operators to return the vectors that build the sparse matrix.
* Modified the m-version of the operators to make them faster.
* Modified matrix assembly in multipatch examples, to make it faster.
* Modified matrix assembly in sp_drchlt_* files.
* sp_nurbs_*_param. Fixed a bug for the functions on the boundary.
* sp_scalar_to_vector_2d. Fixed a bug where the field ndof_dir was missing.

Version geopdes_maxwell-1.0.1

* mp_interface_hcurl_3d.m, fixed a bug in setting the orientation of 
   the functions associated to the edges.

* op_curlu_curlv.cc, fixed a bug where the coefficient was not actually used.

Version geopdes_base-1.0.2 (18/11/2010)

* msh_push_forward_2d.m, msh_push_forward_3d.m. Fixed a bug in the
   computation of the normal vector exterior to the boundary.

* msh_2d_tensor_product.m, msh_3d_tensor_product.m. Added the 
   computation of the normal vector exterior to the boundary.

* sp_scalar_to_vector_3d.m. Fixed a bug where the field ndof_dir 
   was missing.

Version geopdes_base-1.0.1 (17/11/2010)

* inst/space/bsp_2_nrb_1d__.m, inst/space/bsp_2_nrb_2d__.m,
  inst/space/bsp_2_nrb_3d__.m. Fixed a bug where the modified shape
  functions where not returned in the output.

1. INSTALLATION

1.1. OCTAVE

1.1.1 - Install
 
 * Install Octave version 3.2 or higher.  Source code is available from
    http://www.octave.org, binary packages for Mac OSX and Windows are
    available from http://octave.sf.net, Linux binary packages are
    included with all major binary distributions

 * Install the "nurbs" and "integration" packages available on
    http://octave.sf.net: 
   - download the latest release of nurbs-<version>.tar.gz and
      integration-<version>.tar.gz 
   - type, at the octave prompt
      pkg install nurbs-<version>.tar.gz integration-<version>.tar.gz

 * Download the geopdes_base-<version>.tar.gz package from 
    http://geopdes.sf.net/functions/geopdes_base
 
 * Install geopdes_base package by typing at the octave prompt
   pkg install geopdes_base-<version>.tar.gz

 * Load the package:
   pkg load geopdes_base

 * After installing the geopdes_base package, the other packages 
    can be installed in an analogous way.

1.1.3 - Note for windows users

 * The mingw-based binary version of Octave distributed on
    http://octave.sf.net by default installs dynamically loaded binary
    functions in a differnt location as that of script functions. This
    can create problems when installing geopdes_base, it is therefore
    recomended that you change the package installation prefix prior to
    installing geopdes_base, e.g.:

    mkdir c:\octave-packages
    pkg prefix c:\octave-packages c:\octave-packages
    system ("tar xzf geopdes_base-<version>.tar.gz")
    pkg install geopdes_base-<version>

   If this procedure fails for you please let us know, see the section
    "HOW TO CONTRIBUTE" below for the preferred procedure for giving
    feedback. 

1.2. MATLAB

 * Install the "nurbs" package available on http://octave.sf.net:
   - download the latest release of nurbs.tar.gz available on 
      http://octave.sf.net
   - uncompress and untar the file. We recommend you to do this in the
      toolbox folder of Matlab
   - download the C functions for matlab "nurbs_c_files.tar.gz" from
      http://sourceforge.net/projects/geopdes/files
   - uncompress and untar the file in the folder "nurbs/inst"
   - in Matlab, go to the folder "nurbs/inst" and run the script file 'compile'.
      This will compile the files and save the nurbs package to your Matlab path

 * Download the function "grule.m" from
    http://sourceforge.net/projects/geopdes/files
    and save it somewhere in your Matlab path.

 * Install the geopdes_base package:
   - download the latest release of geopdes_base.tar.gz from 
       http://geopdes.sf.net/functions/geopdes_base
   - uncompress and untar the file
   - add the directory "geopdes_base/inst" and its subfolders to
      Matlab's default path (see below)

 * The other packages are installed analogously, but the geopdes_base
    package must be also installed in order to make them work.

 * For convenience, you can set the path using the script "setpath_geopdes.m".
   - download the function from http://sourceforge.net/projects/geopdes/files
   - replace "my_path" by the path of the folder where you saved the packages
   - if necessary, add the path for the function "grule.m" and the nurbs package
   
2. GETTING STARTED

2.1. EXAMPLES
 
 * Each package contains a set of simple examples, that can be run with the 
    command <geopdes_package>_examples. For instance, typing at the Octave prompt

    geopdes_base_examples

    will guide you through a set of menus to choose some simple yet useful 
    examples for solving the Poisson problem.

 * The source code and the data files for these examples can be found in the directory
    <geopdes_package>/inst/examples

 * A collection of simple geometries constructed with NURBS can be found 
    in the following directories

     geopdes_base/inst/examples/geometry_files
     geopdes_multipatch/inst/examples/geometry_files

2.2. HELP

 * A detailed documentation for the software is not likely to be published.
    We suggest you to read the technical report

      C. de Falco, A. Reali, R. Vazquez. 
      GeoPDES: a research tool for Isogeometric Analysis of PDEs, 
      Tech. Report, IMATI-CNR, (2010),

    where we give an explanation of the architecture, the design and
    the main features of the code. This paper can be seen as brief user's guide.

 * The format for the geometry files is explained in the files

     geopdes_base/doc/geo_specs_v07.txt
     geopdes_multipatch/doc/geo_specs_mp_v07.txt

   Notice that multipatch geometries will not run correctly outside the 
    geopdes_multipatch package. At the time of this release, multipatch 
    problems have been also implemented in geopdes_maxwell, geopdes_elasticity,
    and geopdes_fluid for Taylor-Hood elements.

3. HOW TO CONTRIBUTE

3.1 CITING GeoPDEs

 * GeoPDEs has been developed as a part of our research, and is funded by our 
    respective institutions. The best way to support us is by citing

      C. de Falco, A. Reali, R. Vazquez. 
      GeoPDES: a research tool for Isogeometric Analysis of PDEs, 
      Tech. Report, IMATI-CNR, (2010),

    in any paper where GeoPDEs is used to obtain results.

3.2 GIVING FEEDBACK

 * To give us feedback, report bugs, suggest new features or to ask to
    be directly involved in the further development of GeoPDEs, please
    subscribe to the mailing list following the instructions at
    https://lists.sourceforge.net/lists/listinfo/geopdes-users

 *  Before posting a question to the mailing list, please browse the
    mailing list archives at
    https://sourceforge.net/mailarchive/forum.php?forum_name=geopdes-users
    to see whether it has already been answered before.
  
 
