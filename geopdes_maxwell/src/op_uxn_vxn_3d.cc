/* Copyright (C) 2010 Carlo de Falco, Rafael Vazquez

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, write to the Free Software
   Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 02110-1301 USA
*/


#include "geopdes.h"

DEFUN_DLD(op_uxn_vxn_3d, args, nargout,"\n\
OP_UXN_VXN_3D: assemble the matrix M = [m(i,j)], m(i,j) = (mu u_j x n, v_i x n), with n the exterior normal vector. \n\
\n\
  mat = op_uxn_vxn_3d (spu, spv, msh, coeff);\n\
\n\
INPUT:\n\
\n\
 spu:   structure representing the space of trial functions (see sp_bsp_hcurl_3d_phys)\n\
 spv:   structure representing the space of test functions  (see sp_bsp_hcurl_3d_phys)\n\
 msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_3d)\n\
 coeff: physical parameter\n\
\n\
OUTPUT:\n\
\n\
 mat: assembled matrix\n\
")
{
  
  octave_value_list retval;

  geopdes_mesh_normal msh (args(2).map_value ());
  geopdes_space       spu (args(0).map_value (), msh);
  geopdes_space       spv (args(1).map_value (), msh);

  Matrix coeff = args(3).matrix_value ();

  SparseMatrix mat;


  if (!error_state)
    {      
      Array <octave_idx_type>  I (msh.nel () * spv.nsh_max () * spu.nsh_max (), 0.0);
      Array <octave_idx_type> J (msh.nel () * spv.nsh_max () * spu.nsh_max (), 0.0);
      Array <double> V (msh.nel () * spv.nsh_max () * spu.nsh_max (), 0.0);

#pragma omp parallel default (none) shared (msh, spu, spv, I, J, V, coeff)
      {
        double ishp_x_n[3], jshp_x_n[3];
        octave_idx_type counter;

#pragma omp for
        for ( octave_idx_type iel=0; iel < msh.nel (); iel++)
          if (msh.area (iel) > 0.0)
            {
              for ( octave_idx_type idof(0); idof < spv.nsh (iel); idof++) 
                {
                  for ( octave_idx_type jdof(0); jdof < spu.nsh (iel); jdof++)
                    {

                      counter = jdof + spu.nsh (iel) * (idof + spv.nsh (iel) * iel);

                      I(counter) = spv.connectivity (idof, iel)-1;
                      J(counter) = spu.connectivity (jdof, iel)-1;
                      V(counter) = 0.0;

                      for ( octave_idx_type inode(0); inode < msh.nqn (); inode++)
                        {
                          if (msh.weights (inode, iel) > 0.0)
                            {
                              ishp_x_n[0] =
                                spv.shape_functions (1, inode, idof, iel) * msh.normal (2, inode, iel) -
                                spv.shape_functions (2, inode, idof, iel) * msh.normal (1, inode, iel);
                              ishp_x_n[1] =
                                spv.shape_functions (2, inode, idof, iel) * msh.normal (0, inode, iel) -
                                spv.shape_functions (0, inode, idof, iel) * msh.normal (2, inode, iel);
                              ishp_x_n[2] =
                                spv.shape_functions (0, inode, idof, iel) * msh.normal (1, inode, iel) -
                                spv.shape_functions (1, inode, idof, iel) * msh.normal (0, inode, iel);
                              jshp_x_n[0] =
                                spu.shape_functions (1, inode, jdof, iel) * msh.normal (2, inode, iel) -
                                spu.shape_functions (2, inode, jdof, iel) * msh.normal (1, inode, iel);
                              jshp_x_n[1] =
                                spu.shape_functions (2, inode, jdof, iel) * msh.normal (0, inode, iel) -
                                spu.shape_functions (0, inode, jdof, iel) * msh.normal (2, inode, iel);
                              jshp_x_n[2] =
                                spu.shape_functions (0, inode, jdof, iel) * msh.normal (1, inode, iel) -
                                spu.shape_functions (1, inode, jdof, iel) * msh.normal (0, inode, iel);

                              double s = 0.0;
                              for (octave_idx_type icmp(0); icmp < spu.ncomp (); icmp++)
                                s += ishp_x_n[icmp] * jshp_x_n[icmp];

                              V(counter) += msh.jacdet (inode, iel) * msh.weights (inode, iel) * coeff(inode, iel) * s;
                            }  
                        } // end for inode
                      //		   if (idof != jdof) // copy upper triangular part to lower
                      //		      { 
                      //		        I(counter) = J(counter-1);
                      //		        J(counter) = I(counter-1);
                      //		        V(counter) = V(counter-1);
                      //		        counter++;
                      //		      } 
                    } // end for jdof
                } // end for idof
            } else {
#pragma omp critical
            warning_with_id ("geopdes:zero_measure_element", "op_uxn_vxn_3d: element %d has 0 area", iel);
          }  // end for iel, if area > 0
      } // end of parallel section
      mat = SparseMatrix (V, I, J, spv.ndof (), spu.ndof (), true);
      retval (0) = octave_value (mat);
    } // end if !error_state
  return retval;
}

