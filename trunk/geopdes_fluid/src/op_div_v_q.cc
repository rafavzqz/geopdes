/* Copyright (C) 2010, 2011 Carlo de Falco
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

  geopdes_mesh  msh (args(2).scalar_map_value ());
  geopdes_space spv (args(0).scalar_map_value (), msh);
  geopdes_space spq (args(1).scalar_map_value (), msh);

  if (!error_state)
    {

      octave_idx_type iel, inode, idof, jdof, icmp;
      const octave_idx_type nel = msh.nel (), nqn = msh.nqn (), ndof_spq = spq.ndof (), 
        nsh_max_spq = spq.nsh_max (), ndof_spv = spv.ndof (), nsh_max_spv = spv.nsh_max (),
        ncomp_spq = spq.ncomp ();
      
      double shdivv[nsh_max_spv][nqn];
      double shpq[nsh_max_spq][nqn][ncomp_spq];
      double jacdet_weights[nqn];
                              
      octave_idx_type conn_v[nsh_max_spv];
      octave_idx_type conn_q[nsh_max_spq];

      dim_vector dims (nel * nsh_max_spv * nsh_max_spq, 1);

      Array <octave_idx_type> I (dims, 0);
      octave_idx_type* Iptr = I.fortran_vec ();

      Array <octave_idx_type> J (dims, 0);
      octave_idx_type* Jptr = J.fortran_vec ();

      Array <double> V (dims, 0.0);    
      double* Vptr = V.fortran_vec ();

      SparseMatrix mat;

      octave_idx_type counter = 0;
      for ( iel=0; iel < nel; iel++) 
        if (msh.area (iel) > 0.0)
	  {
            for ( inode = 0; inode < nqn; inode++)
              jacdet_weights[inode] = msh.jacdet (inode, iel) * msh.weights (inode, iel);

            // Cache element data to speed-up memory access
            spq.cache_element_shape_functions (iel, (double*)shpq);
            spv.cache_element_shape_function_divs (iel, (double*)shdivv);

            spq.cache_element_connectivity (iel, (octave_idx_type*)conn_q);
            spv.cache_element_connectivity (iel, (octave_idx_type*)conn_v);

            for ( idof = 0; idof <  spq.nsh (iel); idof++) 
              for ( jdof = 0; jdof < spv.nsh (iel); jdof++) 
                {
                  counter = jdof + nsh_max_spv * (idof + nsh_max_spq * iel);
                  
                  Iptr[counter] = conn_q[idof] - 1;
                  Jptr[counter] = conn_v[jdof] - 1;
                  Vptr[counter] = 0.0;

                  for (inode = 0; inode < nqn; inode++)
                    if (msh.weights (inode, iel) > 0.0)
                      Vptr[counter] += jacdet_weights[inode] * shdivv[jdof][inode] * shpq[idof][inode][0];
                }  // end for idof, for jdof
          } 
        else 
          {
            warning_with_id ("geopdes:zero_measure_element", "op_div_v_q: element %d has 0 area (or volume)", iel);
          }  // end for iel, if area > 0

      if (nargout == 1) 
        {
          mat = SparseMatrix (V, I, J, ndof_spq, ndof_spv, true);
          retval(0) = octave_value (mat);
        } 
      else if (nargout == 3)
	{
          for ( icmp = 0; icmp <= counter; icmp++) 
            {
              Iptr[icmp] += 1;
              Jptr[icmp] += 1;
            }
          retval(0) = octave_value (I);
          retval(1) = octave_value (J);
          retval(2) = octave_value (V);
        }
    } // end if !error_state
  return retval;
}
