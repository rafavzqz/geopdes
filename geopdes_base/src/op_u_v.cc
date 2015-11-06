/* Copyright (C) 2009, 2010, 2011 Carlo de Falco
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
      octave_idx_type iel, inode, idof, jdof, icmp;
      octave_idx_type nel = msh.nel ();
      
      const octave_idx_type ncomp = spu.ncomp (), nqn = msh.nqn (),
        ndof_spu = spu.ndof (), nsh_max_spu = spu.nsh_max (), ndof_spv = spv.ndof (), 
        nsh_max_spv = spv.nsh_max ();

      octave_idx_type nsh_u = spu.nsh_max ();
      octave_idx_type nsh_v = spv.nsh_max ();

      double jacdet_weights[nqn];
      double shpv[nsh_v][nqn][ncomp];
      double shpu[nsh_u][nqn][ncomp];
      octave_idx_type conn_v[nsh_v];
      octave_idx_type conn_u[nsh_u];

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

      SparseMatrix mat;

      octave_idx_type counter = 0;
      for (iel=0; iel < nel; iel++)
        if (msh.area (iel) > 0.0)
          {

            for (inode = 0; inode < nqn; inode++)
              jacdet_weights[inode] = msh.jacdet (inode, iel) *
                msh.weights (inode, iel) * coeff (inode, iel);
              

            // Cache element data to speed-up memory access
            spv.cache_element_shape_functions (iel, (double*)shpv);
            spu.cache_element_shape_functions (iel, (double*)shpu);

            spu.cache_element_connectivity (iel, (octave_idx_type*)conn_u);
            spv.cache_element_connectivity (iel, (octave_idx_type*)conn_v);

            for (idof = 0; idof < spv.nsh (iel); idof++) 
              for (jdof = 0; jdof < spu.nsh (iel); jdof++) 
                {
                  counter = jdof + nsh_u * (idof + nsh_v * iel);
                    
                  Iptr[counter] = conn_v[idof] - 1;
                  Jptr[counter] = conn_u[jdof] - 1;
                  Vptr[counter] = 0.0;
                  for (inode = 0; inode < nqn; inode++)
                    if (msh.weights (inode, iel) > 0.0)
                      {
                        double s = 0.0;
                        for (icmp = 0; icmp < ncomp; icmp++)
                          s += shpv[idof][inode][icmp] * 
                            shpu[jdof][inode][icmp];
                          
                        Vptr[counter] += s * jacdet_weights[inode];  
                      }// end for inode, if weight >0  
                } // end for idof, for jdof
          } else {
          warning_with_id ("geopdes:zero_measure_element", 
                           "op_u_v: element %d has 0 area (or volume)", iel);
        }  // end for iel, if area > 0

      if (nargout == 1) 
        {
          mat = SparseMatrix (V, I, J, ndof_spv, ndof_spu, true);
          retval(0) = octave_value (mat);
        } 
      else if (nargout == 3)
	{
          for (icmp = 0; icmp < I.length(); icmp++) 
            {
              Iptr[icmp] += 1;
              Jptr[icmp] += 1;
            }          
          retval(2) = octave_value (V);
          retval(1) = octave_value (J);
          retval(0) = octave_value (I);
        }

    } // end if !error_state
  return retval;
}


