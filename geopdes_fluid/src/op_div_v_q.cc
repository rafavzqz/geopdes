/* Copyright (C) 2010 Carlo de Falco

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

DEFUN_DLD(op_div_v_q, args, nargout,"OP_DIV_V_Q: assemble the matrix B = [b(i,j)], b(i,j) = (q_i, div v_j). \n\
 \n\
   mat = op_div_v_q (spv, spq, msh); \n\
 \n\
 INPUT:  \n\
 \n\
   spv:     structure representing the space of trial functions for the velocity (see sp_bspline_2d_phys)  \n\
   spq:     structure representing the space of test functions for the pressure (see sp_bspline_2d_phys)  \n\
   msh:     structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)  \n\
  \n\
 OUTPUT:  \n\
 \n\
   mat: assembled matrix \n\
")
{
  
  octave_value_list retval;

  geopdes_mesh  msh (args(2).map_value ());
  geopdes_space spv (args(0).map_value (), msh);
  geopdes_space spq (args(1).map_value (), msh);

  if (!error_state)
    {

      Array <octave_idx_type> I (msh.nel () * spv.nsh_max () * spq.nsh_max (), 0);
      Array <octave_idx_type> J (msh.nel () * spv.nsh_max () * spq.nsh_max (), 0);
      Array <double> V (msh.nel () * spv.nsh_max () * spq.nsh_max (), 0.0);      

      SparseMatrix mat;

#pragma omp parallel default (none) shared (msh, spv, spq, I, J, V)
      {

        octave_idx_type counter;
#pragma omp for
      for ( octave_idx_type iel=0; iel < msh.nel (); iel++) 
        if (msh.area (iel) > 0.0)
	  {
	    for ( octave_idx_type idof(0); idof < spq.nsh (iel); idof++) 
	      {
	        for ( octave_idx_type jdof(0); jdof < spv.nsh (iel); jdof++) 
		  {
                    
                    counter = jdof + spv.nsh (iel) * (idof + spq.nsh (iel) * iel);
                    
		    I(counter) = spq.connectivity (idof, iel) - 1;
		    J(counter) = spv.connectivity (jdof, iel) - 1;
		    V(counter) = 0.0;
		    for ( octave_idx_type inode(0); inode < msh.nqn (); inode++)
		      {
		        if (msh.weights (inode, iel) > 0.0)
			  {
			    V(counter) += 
                              msh.jacdet (inode, iel) * msh.weights (inode, iel) *
			      spv.shape_function_divs (inode, jdof, iel) * 
                              spq.shape_functions (0, inode, idof, iel);
			  }  
		      } // end for inode		  
		  } // end for jdof
	      } // end for idof
          } else {
#pragma omp critical
          {warning_with_id ("geopdes:zero_measure_element", "op_div_v_q: element %d has 0 area (or volume)", iel);}
        }  // end for iel, if area > 0
      } // end of parallel region
      mat = SparseMatrix (V, I, J, spq.ndof (), spv.ndof (), true);
      retval(0) = octave_value(mat);
    } // end if !error_state
  return retval;
}
