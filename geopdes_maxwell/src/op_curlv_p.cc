/* Copyright (C) 2009, 2010 Carlo de Falco, Rafael Vazquez

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 2 of the License, or
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

DEFUN_DLD(op_curlv_p, args, nargout,"\n\
OP_CURLV_P: assemble the matrix B = [b(i,j)], b(i,j) = (coeff p_i, curl v_j). \n\
\n\
  mat = op_curlv_p (spu, spv, msh, coeff); \n\
\n\
INPUT:\n\
\n\
  spu:   structure representing the space of vectorial trial functions (see sp_bsp_hcurl_2d_phys) \n\
  spv:   structure representing the space of scalar test functions  (see sp_bsp_l2_2d_phys) \n\
  msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_2d) \n\
  coeff: physical parameter \n\
\n\
OUTPUT:\n\
\n\
  mat: assembled mass matrix")
{
  
  octave_value_list retval;

  geopdes_mesh  msh (args(2).map_value ());
  geopdes_space spv (args(0).map_value (), msh);
  geopdes_space spp (args(1).map_value (), msh);

  Matrix coeff   = args(3).matrix_value ();

  if (!error_state)
    {
      Array <octave_idx_type> I (msh.nel () * spv.nsh_max () * spp.nsh_max (), 0.0);
      Array <octave_idx_type> J (msh.nel () * spv.nsh_max () * spp.nsh_max (), 0.0);
      Array <double> V (msh.nel () * spv.nsh_max () * spp.nsh_max (), 0.0);

      SparseMatrix mat;

#pragma omp parallel default (none) shared (msh, spp, spv, I, J, V)
      {      
        octave_idx_type counter;

#pragma omp for
        for ( octave_idx_type iel=0; iel < msh.nel (); iel++) 
          if (msh.area (iel) > 0.0)
            {
              for ( octave_idx_type idof(0); idof < spp.nsh (iel); idof++) 
                {
                  for (octave_idx_type jdof(0); jdof < spv.nsh (iel); jdof++)
                    {
                      counter = jdof + spv.nsh (iel) * (idof + spp.nsh (iel) * iel);
                      I(counter) = spp.connectivity (idof,iel) - 1;
                      J(counter) = spv.connectivity (jdof,iel) - 1;
                      V(counter) = 0.0;
                      for ( octave_idx_type inode(0); inode < msh.nqn (); inode++)
                        {
                          if (msh.weights (inode, iel) > 0.0)
                            {			  
                              V(counter) += 
                                msh.jacdet (inode, iel) * msh.weights (inode, iel) * 
                                (spp.shape_functions (0, inode, idof, iel) * spv.shape_function_curls (inode, jdof, iel));
                            }  
                        } // end for inode		  
                      counter++;
                    } // end for jdof
                } // end for idof
            } else {
            warning_with_id ("geopdes:zero_measure_element", "op_curlv_p: element %d has 0 area", iel);
          }  // end for iel, if area > 0
      } // end parallel section
      mat = SparseMatrix(V, I, J, spp.ndof (), spv.ndof (), true);
      retval(0) = octave_value (mat);
    } // end if !error_state
  return retval;
}
