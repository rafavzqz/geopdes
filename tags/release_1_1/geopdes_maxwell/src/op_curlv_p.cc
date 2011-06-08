/* Copyright (C) 2009, 2010 Carlo de Falco, Rafael Vazquez
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

DEFUN_DLD(op_curlv_p, args, nargout,"\n\
OP_CURLV_P: assemble the matrix B = [b(i,j)], b(i,j) = (coeff p_i, curl v_j). \n\
\n\
  mat = op_curlv_p (spv, spp, msh, coeff); \n\
\n\
INPUT:\n\
\n\
  spv:   structure representing the space of vectorial trial functions (see sp_bsp_hcurl_2d_phys) \n\
  spp:   structure representing the space of scalar test functions  (see sp_bsp_l2_2d_phys) \n\
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
      const octave_idx_type nel = msh.nel (), ndir = msh.ndir (), ncomp = spv.ncomp (), nqn = msh.nqn (), ndof_spp = spp.ndof (), nsh_max_spp = spp.nsh_max (), ndof_spv = spv.ndof (), nsh_max_spv = spv.nsh_max ();

      dim_vector dims (nel * nsh_max_spv * nsh_max_spp, 1);
      Array <octave_idx_type> I (dims, 0);
      Array <octave_idx_type> J (dims, 0);
      Array <double> V (dims, 0.0);

      SparseMatrix mat;

      octave_idx_type counter = 0, iel, inode, idof, jdof, icmp;

#pragma omp parallel default (none) shared (msh, spp, spv, I, J, V)
      {      
#pragma omp for
        for ( iel=0; iel < nel; iel++) 
          if (msh.area (iel) > 0.0)
            {
              const octave_idx_type nsh_p = spp.nsh (iel);
              const octave_idx_type nsh_v = spv.nsh (iel);
              double jacdet_weights[nqn];

              for ( inode = 0; inode < nqn; inode++)
                {
                  jacdet_weights[inode] = msh.jacdet (inode, iel) *
                    msh.weights (inode, iel) * coeff (inode, iel);
                }

              double shcv[nsh_v][nqn];
              double shpp[nsh_p][nqn];
              int conn_v[nsh_v];
              int conn_p[nsh_p];

              for ( idof = 0; idof < nsh_p; idof++) 
		{
                  for ( inode = 0; inode < nqn; inode++)
                    {
                      shpp[idof][inode] = spp.shape_functions (0, inode, idof, iel);
                    }
		  conn_p[idof] = spp.connectivity (idof, iel) - 1;
	        }

              for ( jdof = 0; jdof < nsh_v; jdof++) 
		{
                  for ( inode = 0; inode < nqn; inode++)
                    {
                      shcv[jdof][inode] = spv.shape_function_curls (inode, jdof, iel);
                    }
		  conn_v[jdof] = spv.connectivity (jdof, iel) - 1;
	        }


              for ( idof = 0; idof < nsh_p; idof++) 
                {
                  for ( jdof = 0; jdof < nsh_v; jdof++) 
                    {
                      counter = jdof + nsh_v * (idof + nsh_p * iel);

                      I(counter) = conn_p[idof];
                      J(counter) = conn_v[jdof];
                      V(counter) = 0.0;
                      for ( inode = 0; inode < nqn; inode++)
                        {
                          if (msh.weights (inode, iel) > 0.0)
                            {			  
                              V(counter) += jacdet_weights[inode] * 
                                shpp[idof][inode] * shcv[jdof][inode];
                            }  
                        } // end for inode		  
                      counter++;
                    } // end for jdof
                } // end for idof
            } else {
#pragma omp critical
            {warning_with_id ("geopdes:zero_measure_element", "op_curlv_p: element %d has 0 area", iel);}
          }  // end for iel, if area > 0
      } // end parallel section

      if (nargout == 1) 
        {
          mat = SparseMatrix (V, I, J, ndof_spp, ndof_spv, true);
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
