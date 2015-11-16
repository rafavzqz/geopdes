/* Copyright (C) 2009, 2010, 2011 Carlo de Falco
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

DEFUN_DLD(op_curlu_curlv_3d, args, nargout,"\n\
OP_CURLU_CURLV_3D: assemble the matrix A = [a(i,j)], a(i,j) = (coeff curl u_j, curl v_i), with a vector-valued curl. \n\
\n\
  mat = op_curlu_curlv_3D (spu, spv, msh, epsilon); \n\
\n\
INPUT: \n\
\n\
  spu:   structure representing the space of trial functions (see sp_vector/sp_evaluate_col) \n\
  spv:   structure representing the space of test functions  (see sp_vector/sp_evaluate_col) \n\
  msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col) \n\
  coeff: physical parameter \n\
\n\
OUTPUT: \n\
\n\
  mat: assembled matrix \n\
")
{
  
  octave_value_list retval;

  geopdes_mesh  msh (args(2).scalar_map_value ());
  geopdes_space spu (args(0).scalar_map_value (), msh);
  geopdes_space spv (args(1).scalar_map_value (), msh);
  Matrix coeff   = args(3).matrix_value();

  if (!error_state)
    {

      const octave_idx_type nel = msh.nel (), nqn = msh.nqn (), ndof_spu = spu.ndof (), nsh_max_spu = spu.nsh_max (), ndof_spv = spv.ndof (), nsh_max_spv = spv.nsh_max ();

      dim_vector dims (nel * nsh_max_spv * nsh_max_spu, 1);
      Array <octave_idx_type> I (dims, 0);
      octave_idx_type* Iptr = I.fortran_vec ();

      Array <octave_idx_type> J (dims, 0);
      octave_idx_type* Jptr = J.fortran_vec ();

      Array <double> V (dims, 0.0);
      double* Vptr = V.fortran_vec ();

      SparseMatrix mat;

      octave_idx_type counter = 0, iel, inode, idof, jdof, icmp;

      for (iel=0; iel < nel; iel++) 
        if (msh.area (iel) > 0.0)
	  {
            const octave_idx_type nsh_u = spu.nsh (iel);
            const octave_idx_type nsh_v = spv.nsh (iel);
            double jacdet_weights[nqn];

            for (inode = 0; inode < nqn; inode++)
              jacdet_weights[inode] = msh.jacdet (inode, iel) *
                msh.weights (inode, iel) * coeff (inode, iel);

            double shcv[nsh_v][nqn][3];
            double shcu[nsh_u][nqn][3];
            octave_idx_type conn_v[nsh_max_spv];
            octave_idx_type conn_u[nsh_max_spu];

            for (idof = 0; idof < nsh_v; idof++) 
              for (inode = 0; inode < nqn; inode++)
                for (icmp = 0; icmp < 3; icmp++)
                  shcv[idof][inode][icmp] = spv.shape_function_curls (icmp, inode, idof, iel);
            
            for (jdof = 0; jdof < nsh_u; jdof++) 
              for (inode = 0; inode < nqn; inode++)
                for (icmp = 0; icmp < 3; icmp++)
                  shcu[jdof][inode][icmp] = spu.shape_function_curls (icmp, inode, jdof, iel);
            
            spu.cache_element_connectivity (iel, (octave_idx_type*)conn_u);
            spv.cache_element_connectivity (iel, (octave_idx_type*)conn_v);

            for (idof = 0; idof < nsh_v; idof++) 
              for (jdof = 0; jdof < nsh_u; jdof++) 
                {
                  counter = jdof + nsh_u * (idof + nsh_v * iel);
                  
                  Iptr[counter] = conn_v[idof] - 1;
                  Jptr[counter] = conn_u[jdof] - 1;
                  Vptr[counter] = 0.0;
                  for (inode = 0; inode < nqn; inode++)
                    if (msh.weights (inode, iel) > 0.0)
                      {			  
                        double s = 0.0;
                        for (icmp = 0; icmp < 3; icmp++)
                          s += shcv[idof][inode][icmp] * 
                            shcu[jdof][inode][icmp];
                        
                          Vptr[counter] += jacdet_weights[inode] * s;
                      }  // end for inode, if 
                } // end for idof, for jdof
          } 
        else
          {
            warning_with_id ("geopdes:zero_measure_element", "op_curlu_curlv_3d: element %d has 0 volume", iel);
          }  // end for iel, if area > 0

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
