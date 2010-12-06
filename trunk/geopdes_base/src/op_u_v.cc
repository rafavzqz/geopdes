/* Copyright (C) 2009 Carlo de Falco

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

DEFUN_DLD(op_u_v, args, nargout,"\n\
OP_U_V: assemble the mass matrix M = [m(i,j)], m(i,j) = (mu u_j, v_i).\n\
\n\
  mat = op_u_v (spu, spv, msh, coeff); \n\
\n\
INPUT: \n\
\n\
  spu:   structure representing the space of trial functions (see sp_bspline_2d_phys) \n\
  spv:   structure representing the space of test functions  (see sp_bspline_2d_phys) \n\
  msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_2d) \n\
  coeff: reaction coefficient \n\
\n\
OUTPUT: \n\
\n\
  mat: assembled mass matrix \n\
")
{
  
  octave_value_list retval;

  geopdes_mesh  msh (args(2).map_value ());
  geopdes_space spu (args(0).map_value (), msh);
  geopdes_space spv (args(1).map_value (), msh);

  Matrix coeff = args(3).matrix_value ();

  if (!error_state)
    {      

      ColumnVector I (msh.nel () * spv.nsh_max () * spu.nsh_max (), 0.0);
      ColumnVector J (msh.nel () * spv.nsh_max () * spu.nsh_max (), 0.0);
      ColumnVector V (msh.nel () * spv.nsh_max () * spu.nsh_max (), 0.0);

      SparseMatrix mat;

#pragma omp parallel default (none) shared (msh, spu, spv, I, J, V, coeff)
      {
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
                              double s = 0.0;
                              for (octave_idx_type icmp(0); icmp < spu.ncomp (); icmp++)
                                s += spv.shape_functions (icmp, inode, idof, iel) *
                                  spu.shape_functions (icmp, inode, jdof, iel);
                              
                              V(counter) += msh.jacdet (inode, iel) * 
                                msh.weights (inode, iel) * 
                                coeff(inode, iel) * s;			  

                            }  
                        } // end for inode
                    } // end for jdof
                } // end for idof
            } else {
            warning_with_id ("geopdes:zero_measure_element", "op_u_v: element %d has 0 area (or volume)", iel);
          }  // end for iel, if area > 0
      } // end of openmp parallel section

      mat = SparseMatrix (V, I, J, spv.ndof (), spu.ndof (), true);
      retval (0) = octave_value (mat);

    } // end if !error_state
  return retval;
}


