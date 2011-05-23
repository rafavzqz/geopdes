/* Copyright (C) 2009, 2010 Carlo de Falco
   Copyright (C) 2011 Rafael Vazquez

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

DEFUN_DLD(op_curlu_curlv_2d, args, nargout,"\n\
OP_CURLU_CURLV_2D: assemble the matrix A = [a(i,j)], a(i,j) = (coeff curl u_j, curl v_i), with a scalar-valued curl. \n\
\n\
  mat = op_curlu_curlv_2d (spu, spv, msh, epsilon); \n\
\n\
INPUT: \n\
\n\
  spu:   structure representing the space of trial functions (see sp_bsp_Hcurl_2d_phys) \n\
  spv:   structure representing the space of test functions  (see sp_bsp_Hcurl_2d_phys) \n\
  msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_2d) \n\
  coeff: physical parameter \n\
\n\
OUTPUT: \n\
\n\
  mat: assembled matrix \n\
")
{
  
  octave_value_list retval;
  
  geopdes_mesh  msh (args(2).map_value ());
  geopdes_space spu (args(0).map_value (), msh);
  geopdes_space spv (args(1).map_value (), msh);
  
  Matrix coeff   = args(3).matrix_value();

  if (!error_state)
    {
      const octave_idx_type nel = msh.nel (), nqn = msh.nqn (), ndof_spu = spu.ndof (), nsh_max_spu = spu.nsh_max (), ndof_spv = spv.ndof (), nsh_max_spv = spv.nsh_max ();

      dim_vector dims (nel * nsh_max_spv * nsh_max_spu, 1);
      Array <octave_idx_type> I (dims, 0);
      Array <octave_idx_type> J (dims, 0);
      Array <double> V (dims, 0.0);
      
      SparseMatrix mat;

      octave_idx_type counter = 0, iel, inode, idof, jdof;

#pragma omp parallel default (none) shared (msh, spu, spv, I, J, V, coeff)
      {      
#pragma omp for
      for ( iel=0; iel < nel; iel++) 
        if (msh.area (iel) > 0.0)
	  {
            const octave_idx_type nsh_u = spu.nsh (iel);
            const octave_idx_type nsh_v = spv.nsh (iel);
            double jacdet_weights[nqn];

            for ( inode = 0; inode < nqn; inode++)
              {
                jacdet_weights[inode] = msh.jacdet (inode, iel) *
                  msh.weights (inode, iel) * coeff (inode, iel);
              }

            double shcv[nsh_v][nqn];
            double shcu[nsh_u][nqn];
            int conn_v[nsh_v];
            int conn_u[nsh_u];

              for ( idof = 0; idof < nsh_v; idof++) 
		{
                  for ( inode = 0; inode < nqn; inode++)
                    {
                      shcv[idof][inode] = spv.shape_function_curls (inode, idof, iel);
                    }
		  conn_v[idof] = spv.connectivity (idof, iel) - 1;
	        }

              for ( jdof = 0; jdof < nsh_u; jdof++) 
		{
                  for ( inode = 0; inode < nqn; inode++)
                    {
                      shcu[jdof][inode] = spu.shape_function_curls (inode, jdof, iel);
                    }
		  conn_u[jdof] = spu.connectivity (jdof, iel) - 1;
	        }

            for ( idof = 0; idof < nsh_v; idof++) 
	      {
                for ( jdof = 0; jdof < nsh_u; jdof++) 
		  {
                    counter = jdof + nsh_u * (idof + nsh_v * iel);
                    
                    I(counter) = conn_v[idof];
                    J(counter) = conn_u[jdof];
		    V(counter) = 0.0;
                    for ( inode = 0; inode < nqn; inode++)
		      {
		        if (msh.weights (inode, iel) > 0.0)
			  {			  
			    V(counter) += jacdet_weights[inode] * 
			      shcv[idof][inode] * shcu[jdof][inode];
			  }  
		      } // end for inode		  
//		    if (idof != jdof) // copy upper triangular part to lower
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
          {warning_with_id ("geopdes:zero_measure_element", "op_curlu_curlv_2d: element %d has 0 area", iel);}
        }  // end for iel, if area > 0
      } // end of parallel section
      mat = SparseMatrix (V, I, J, ndof_spv, ndof_spu, true);
      retval(0) = octave_value (mat);
    } // end if !error_state
  return retval;
}
