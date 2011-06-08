/* Copyright (C) 2010 Carlo de Falco
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

      const octave_idx_type nel = msh.nel (), nqn = msh.nqn (), ndof_spq = spq.ndof (), nsh_max_spq = spq.nsh_max (), ndof_spv = spv.ndof (), nsh_max_spv = spv.nsh_max ();

      dim_vector dims (nel * nsh_max_spv * nsh_max_spq, 1);
      Array <octave_idx_type> I (dims, 0);
      Array <octave_idx_type> J (dims, 0);
      Array <double> V (dims, 0.0);    

      SparseMatrix mat;

      octave_idx_type counter = 0, iel, inode, idof, jdof, icmp;

#pragma omp parallel default (none) shared (msh, spv, spq, I, J, V)
      {

#pragma omp for
      for ( iel=0; iel < nel; iel++) 
        if (msh.area (iel) > 0.0)
	  {
            const octave_idx_type nsh_q = spq.nsh (iel);
            const octave_idx_type nsh_v = spv.nsh (iel);
            double jacdet_weights[nqn];

            for ( inode = 0; inode < nqn; inode++)
              {
                jacdet_weights[inode] = msh.jacdet (inode, iel) *
                  msh.weights (inode, iel);
              }

            double shdivv[nsh_v][nqn];
            double shpq[nsh_q][nqn];
            int conn_v[nsh_v];
            int conn_q[nsh_q];

              for ( jdof = 0; jdof < nsh_v; jdof++)
		{
                  for ( inode = 0; inode < nqn; inode++)
                    {
                      shdivv[jdof][inode] = spv.shape_function_divs (inode, jdof, iel);
                    }
		  conn_v[jdof] = spv.connectivity (jdof, iel) - 1;
	        }

              for ( idof = 0; idof < nsh_q; idof++)
		{
                  for ( inode = 0; inode < nqn; inode++)
                    {
                      shpq[idof][inode] = spq.shape_functions (0, inode, idof, iel);
                    }
		  conn_q[idof] = spq.connectivity (idof, iel) - 1;
	        }


            for ( idof = 0; idof < nsh_q; idof++) 
	      {
                for ( jdof = 0; jdof < nsh_v; jdof++) 
		  {
                    counter = jdof + nsh_v * (idof + nsh_q * iel);
                    
                    I(counter) = conn_q[idof];
                    J(counter) = conn_v[jdof];
		    V(counter) = 0.0;
                    for ( inode = 0; inode < nqn; inode++)
		      {
		        if (msh.weights (inode, iel) > 0.0)
			  {
			    V(counter) += jacdet_weights[inode] *
                              shdivv[jdof][inode] * shpq[idof][inode];
			  }  
		      } // end for inode		  
		  } // end for jdof
	      } // end for idof
          } else {
#pragma omp critical
          {warning_with_id ("geopdes:zero_measure_element", "op_div_v_q: element %d has 0 area (or volume)", iel);}
        }  // end for iel, if area > 0
      } // end of parallel region

      if (nargout == 1) 
        {
          mat = SparseMatrix (V, I, J, ndof_spq, ndof_spv, true);
          retval(0) = octave_value (mat);
        } 
      else if (nargout == 3)
	{
          for ( icmp = 0; icmp <= counter; icmp++) 
            {
              I(icmp)++;
              J(icmp)++;
            }
          retval(0) = octave_value (I);
          retval(1) = octave_value (J);
          retval(2) = octave_value (V);
        }
    } // end if !error_state
  return retval;
}
