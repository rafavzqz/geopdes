/* Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
   Copyright (C) 2011 Rafael Vazquez

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
      const octave_idx_type nel = msh.nel (), ncomp = spu.ncomp (), nqn = msh.nqn (), ndof_spu = spu.ndof (), nsh_max_spu = spu.nsh_max (), ndof_spv = spv.ndof (), nsh_max_spv = spv.nsh_max ();

      dim_vector dims (nel * nsh_max_spv * nsh_max_spu, 1);
      Array <octave_idx_type> I (dims, 0);
      Array <octave_idx_type> J (dims, 0);
      Array <double> V (dims, 0.0);

#pragma omp parallel default (none) shared (msh, spu, spv, I, J, V, coeff)
      {
        octave_idx_type counter = 0, iel, inode, idof, jdof, icmp;
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

              double shpv_x_n[nsh_v][nqn][ncomp];
              double shpu_x_n[nsh_u][nqn][ncomp];
              int conn_v[nsh_v];
              int conn_u[nsh_u];

              for ( idof = 0; idof < nsh_v; idof++) 
		{
                  for ( inode = 0; inode < nqn; inode++)
                    {
                      for ( icmp = 0; icmp < ncomp; icmp++)
                        {
                          shpv_x_n[idof][inode][0] =
                                spv.shape_functions (1, inode, idof, iel) * 
                                msh.normal (2, inode, iel) -
                                spv.shape_functions (2, inode, idof, iel) * 
                                msh.normal (1, inode, iel);
                          shpv_x_n[idof][inode][1] =
                                spv.shape_functions (2, inode, idof, iel) * 
                                msh.normal (0, inode, iel) -
                                spv.shape_functions (0, inode, idof, iel) * 
                                msh.normal (2, inode, iel);
                          shpv_x_n[idof][inode][2] =
                                spv.shape_functions (0, inode, idof, iel) * 
                                msh.normal (1, inode, iel) -
                                spv.shape_functions (1, inode, idof, iel) * 
                                msh.normal (0, inode, iel);
                        }
                    }
		  conn_v[idof] = spv.connectivity (idof, iel) - 1;
	        }

              for ( jdof = 0; jdof < nsh_u; jdof++) 
		{
                  for ( inode = 0; inode < nqn; inode++)
                    {
                      for ( icmp = 0; icmp < ncomp; icmp++)
                        {
                          shpu_x_n[jdof][inode][0] =
                                spu.shape_functions (1, inode, jdof, iel) * 
                                msh.normal (2, inode, iel) -
                                spu.shape_functions (2, inode, jdof, iel) * 
                                msh.normal (1, inode, iel);
                          shpu_x_n[jdof][inode][1] =
                                spu.shape_functions (2, inode, jdof, iel) * 
                                msh.normal (0, inode, iel) -
                                spu.shape_functions (0, inode, jdof, iel) * 
                                msh.normal (2, inode, iel);
                          shpu_x_n[jdof][inode][2] =
                                spu.shape_functions (0, inode, jdof, iel) * 
                                msh.normal (1, inode, iel) -
                                spu.shape_functions (1, inode, jdof, iel) * 
                                msh.normal (0, inode, iel);
                        }
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
                              double s = 0.0;
                              for ( icmp = 0; icmp < ncomp; icmp++)
                                s += shpv_x_n[idof][inode][icmp] * 
                                  shpu_x_n[jdof][inode][icmp];

                              V(counter) += s * jacdet_weights[inode];
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
            {warning_with_id ("geopdes:zero_measure_element", "op_uxn_vxn_3d: element %d has 0 area", iel);}
          }  // end for iel, if area > 0
      } // end of parallel section
      mat = SparseMatrix (V, I, J, ndof_spv, ndof_spu, true);
      retval (0) = octave_value (mat);
    } // end if !error_state
  return retval;
}

