/* Copyright (C) 2010, 2011 Carlo de Falco, Rafael Vazquez

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
      const octave_idx_type nel = msh.nel (), ncomp = spu.ncomp (), nqn = msh.nqn (), 
        ndof_spu = spu.ndof (), nsh_max_spu = spu.nsh_max (), ndof_spv = spv.ndof (), nsh_max_spv = spv.nsh_max ();

#if OCTAVE_API_VERSION_NUMBER>37
      dim_vector dims (nel * nsh_max_spv * nsh_max_spu, 1);
      Array <octave_idx_type> I (dims, 0);
      octave_idx_type* Iptr = I.fortran_vec ();

      Array <octave_idx_type> J (dims, 0);
      octave_idx_type* Jptr = J.fortran_vec ();

      Array <double> V (dims, 0.0);    
      double* Vptr = V.fortran_vec ();
#else
      ColumnVector I (nel * nsh_max_spv * nsh_max_spu, 0);
      double* Iptr = I.fortran_vec ();

      ColumnVector J (nel * nsh_max_spv * nsh_max_spu, 0);
      double* Jptr = J.fortran_vec ();

      ColumnVector V (nel * nsh_max_spv * nsh_max_spu, 0.0);
      double* Vptr = V.fortran_vec ();
#endif

      octave_idx_type counter = 0, iel, inode, idof, jdof, icmp;

      {
        for (iel=0; iel < nel; iel++) 
          if (msh.area (iel) > 0.0)
            {
              const octave_idx_type nsh_u = spu.nsh (iel);
              const octave_idx_type nsh_v = spv.nsh (iel);
              double jacdet_weights[nqn];

              for (inode = 0; inode < nqn; inode++)
                {
                  jacdet_weights[inode] = msh.jacdet (inode, iel) *
                    msh.weights (inode, iel) * coeff (inode, iel);
                }

              double shpv_x_n[nsh_v][nqn][ncomp];
              double shpu_x_n[nsh_u][nqn][ncomp];
              octave_idx_type conn_v[nsh_v];
              octave_idx_type conn_u[nsh_u];

              // Cache element data to speed-up memory access
              spu.cache_element_connectivity (iel, (octave_idx_type*)conn_u);
              spv.cache_element_connectivity (iel, (octave_idx_type*)conn_v);

              for (idof = 0; idof < nsh_v; idof++) 
                for (inode = 0; inode < nqn; inode++)
                  for (icmp = 0; icmp < ncomp; icmp++)
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

              for (jdof = 0; jdof < nsh_u; jdof++) 
                for (inode = 0; inode < nqn; inode++)
                  for (icmp = 0; icmp < ncomp; icmp++)
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
              
              for (idof = 0; idof < nsh_v; idof++) 
                {
                  for (jdof = 0; jdof < nsh_u; jdof++) 
                    {
                      counter = jdof + nsh_u * (idof + nsh_v * iel);

                      Iptr[counter] = conn_v[idof] - 1;
                      Jptr[counter] = conn_u[jdof] - 1;
                      Vptr[counter] = 0.0;
                      for (inode = 0; inode < nqn; inode++)
                        {
                          if (msh.weights (inode, iel) > 0.0)
                            {
                              double s = 0.0;
                              for (icmp = 0; icmp < ncomp; icmp++)
                                s += shpv_x_n[idof][inode][icmp] * 
                                  shpu_x_n[jdof][inode][icmp];

                              Vptr[counter] += s * jacdet_weights[inode];
                            }  
                        } // end for inode
                    } // end for jdof
                } // end for idof
            } else {
            {warning_with_id ("geopdes:zero_measure_element", "op_uxn_vxn_3d: element %d has 0 area", iel);}
          }  // end for iel, if area > 0
      } // end of parallel section

      if (nargout == 1) 
        {
          mat = SparseMatrix (V, I, J, ndof_spv, ndof_spu, true);
          retval(0) = octave_value (mat);
        } 
      else if (nargout == 3)
	{
          for (icmp = 0; icmp <= counter; icmp++) 
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

