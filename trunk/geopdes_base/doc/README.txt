CONTENTS:
1: CONTENTS AND INSTALLATION
  1.1: OCTAVE
    1.1.1 - Install
    1.1.2 - Note for Windows users
  1.2: MATLAB
2: GETTING STARTED
  2.1: EXAMPLES
  2.2: HELP
3: HOW TO CONTRIBUTE
  3.1: CITING GeoPDEs
  3.2: PROVIDING FEEDBACK
4: CHANGELOG

-----------------------------------------------------

1. CONTENTS AND INSTALLATION

  Download and uncompress the file GeoPDEs_full.tar.gz. This file contains 
   all that is necessary to install GeoPDEs:
   - The NURBS toolbox: nurbs-1.3.6.tar.gz
   - The GeoPDEs packages:
      geopdes_base-2.0.3.tar.gz, geopdes_elasticity-2.0.2.tar.gz, 
      geopdes_fluid-2.0.2.tar.gz, geopdes_maxwell-2.0.2.tar.gz, 
      geopdes_multipatch-2.0.3.tar.gz, geopdes_tsplines-2.0.3.tar.gz
   - The technical report: GeoPDES_report.pdf
   - This README file.
   - For Matlab users, the mex files of the NURBS toolbox: nurbs_mex_files.tar.gz
   - For Matlab users, the file geopdes_setpath.m
   
    
1.1. OCTAVE

1.1.1 - Install
 
 * Be sure to have installed Octave version 3.2 or higher. 
    Source code is available from http://www.octave.org, 
    binary packages for Mac OSX and Windows are available
    from http://octave.sf.net, Linux binary packages are
    included with all major binary distributions.

 * Install the "nurbs" package by typing at the octave prompt
    pkg install nurbs-<version>.tar.gz
    pkg load nurbs

 * Install geopdes_base package by typing at the octave prompt
    pkg install geopdes_base-<version>.tar.gz
    pkg load geopdes_base

 * After installing the geopdes_base package, the other packages 
    can be installed in an analogous way.

1.1.2 - Note for Windows users

   The mingw-based binary version of Octave distributed on
    http://octave.sf.net by default installs dynamically loaded binary
    functions in a different location as that of script functions. This
    can create problems when installing geopdes_base. The procedure to 
    install the package in Windows is as follows: 

 * Install octave and the nurbs package, like in
    steps 1 to 2 of Section 1.1.1.

 * Be sure to put the file geopdes_base-<version>.tar.gz to a folder 
    without blanks in the name, for instance, c:\tmp

 * In Octave, change the package installation prefix and install geopdes_base:
     cd c:\tmp
     mkdir c:\octave-packages
     pkg prefix c:\octave-packages c:\octave-packages
     system ("tar xzf geopdes_base-<version>.tar.gz")
     pkg install geopdes_base

   If this procedure fails for you please let us know. See the section
    "HOW TO CONTRIBUTE" below for the preferred procedure for giving
    feedback.

1.2. MATLAB

 * Install the "nurbs" package available on http://octave.sf.net:
   - uncompress and untar the file nurbs-<version>.tar.gz. 

 * Install the mex-files for the "nurbs" package (OPTIONAL):
   - uncompress and untar the file nurbs_mex_files.tar.gz in the folder 
      "nurbs/inst" of the nurbs package
   - in Matlab, go to the folder "nurbs/inst" and run the script file 'compile'.
      This will compile the files and save the nurbs package to your Matlab path

 * Install the geopdes_base package:
   - uncompress and untar the file geopdes_base-<version>.tar.gz
   - add the directory "geopdes_base/inst" and its subfolders to
      Matlab's default path (see below)

 * The other packages are installed analogously, but the geopdes_base
    package must be also installed in order to make them work.

 * For convenience, you can set the path using the script "geopdes_setpath.m".
   - replace "my_path" by the path of the folder where you saved the packages
   - run the script, typing "geopdes_setpath" in Matlab command window
   
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
    We suggest you to read the paper

      C. de Falco, A. Reali, R. Vazquez. 
      GeoPDES: a research tool for Isogeometric Analysis of PDEs, 
      Advances in Engineering Software, 42(12), 2011
      doi:10.1016/j.advengsoft.2011.06.10

    also available in preprint form as

      Tech. Report, IMATI-CNR, (2010),

    In the paper we give an explanation of the architecture, the design and
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

 * GeoPDEs has been developed as a part of our research, and is funded by 
    our respective institutions. The best way to support us is by citing

      C. de Falco, A. Reali, R. Vazquez. 
      GeoPDES: a research tool for Isogeometric Analysis of PDEs, 
      Advances in Engineering Software, 42(12), 2011
      doi:10.1016/j.advengsoft.2011.06.010

    in any paper where GeoPDEs is used to obtain numerical results.

3.2 GIVING FEEDBACK

 * To give us feedback, report bugs, suggest new features or to ask to
    be directly involved in the further development of GeoPDEs, please
    subscribe to the mailing list following the instructions at
    https://lists.sourceforge.net/lists/listinfo/geopdes-users

 * Before posting a question to the mailing list, please browse the
    mailing list archives at
    https://sourceforge.net/mailarchive/forum.php?forum_name=geopdes-users
    to see whether it has already been answered before.
  
4: CHANGELOG

Below is a list of the main changes introduced with each new release

Version geopdes_***-2.0.4
* geopdes_base, added functions for the advection-diffusion problem with SUPG stabilization: op_mat_stab_SUPG, op_mat_stab_SUPG_tp, op_rhs_stab_SUPG, op_rhs_stab_SUPG_tp, solve_adv_diff_2d, ex_advection_diffusion_square.
* geopdes_base, added functions for the bilaplacian: op_laplaceu_laplacev, op_laplaceu_laplacev_tp, solve_bilaplace_2d_iso.
* geopdes_elasticity, added functions for the bilaplacian: op_gradgradu_gradgradv, op_gradgradu_gradgradv_tp, solve_bilaplace_gradgrad_2d_iso, ex_kirchhoff_*.
* geopdes_fluid, added functions for the Nitsche's method: sp_weak_dirichlet_bc, op_udotn_vdotn, op_fdotn_vdotn, op_gradv_n_u, op_gradv_n_f.
* solve_stokes_2d, sp_bspline_fluid_2d: modified to use Nitsche's method.
* msh_2d/msh_eval_boundary_side: added computation of the normal characteristic length.
* sp_eval_stress: modified to work also with the Piola transformation.

Version geopdes_***-2.0.3

* In geopdes.h, "quad_nodes" replaced by "geo_map_jac" to check the dimension.
* Fixed bug in the functions mp_interface_#d

Version geopdes_***-2.0.2

* sp_eval and sp_to_vtk: added the possibility to plot curl and gradient
* Added the function(s) sp_eval_grad(div,curl)_msh
* Added the old version of the operators, at the end of the new ones,
    for didactic purposes
* Linear elasticity: added pressure and symmetric conditions in 3D
* Linear elasticity: changed the name of the variables for Lame parameters
* @msh_3d/msh_precompute: fixed bug in the loop

Version geopdes_***-2.0.1

* Modified the oct-operators to be compatible with Octave 3.2
* Fixed bug in the function names sp_vector_#d_curl_transform

Version geopdes_***-2.0.0

* Added the classes msh_2d and msh_3d, with their methods.
* Added the classes sp_bspline_2d, sp_bspline_3d, sp_nurbs_2d, sp_nurbs_3d,
    sp_vector_2d, sp_vector_3d, sp_vector_2d_piola_transform, 
    sp_vector_2d_curl_transform, sp_vector_3d_curl_transform,
    sp_bspline_2d_3forms, with their methods.
* Removed many functions that became unnecessary with the classes.
* Functions sp_eval and sp_to_vtk work now for 2D and 3D geometries.
* Added the functions sp_eval_msh, sp_eval_div_msh.
* Added the tensor product version of the operators (like op_u_v_tp.m).
    They only work for classes.
* Modified the m-version of the operators, to make them even faster.

Version geopdes_***-1.1.0

* Added the function grule to the "utils" folder.
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

Version geopdes_maxwell-1.0.1 (not released)

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
