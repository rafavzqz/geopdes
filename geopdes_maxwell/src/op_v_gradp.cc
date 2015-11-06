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

DEFUN_DLD(op_v_gradp, args, nargout,"\n\
OP_V_GRADP: assemble the matrix B = [b(i,j)], b(i,j) = (epsilon grad p_i, v_j). \n\
\n\
   mat = op_v_gradp (spv, spp, msh, epsilon);\n\
\n\
 INPUT:\n\
\n\
   spv:     structure representing the space of vectorial trial functions  (see sp_bsp_hcurl_2d_phys). \n\
   spp:     structure representing the space of scalar test functions (see sp_bspline_2d_phys). \n\
   msh:     structure containing the domain partition and the quadrature rule (see msh_push_forward_2d). \n\
   epsilon: physical parameter. \n\
\n\
 OUTPUT: \n\
\n\
   mat: assembled matrix. \n\
")
{
  
  octave_value_list retval;

  geopdes_mesh  msh (args(2).map_value ());
  geopdes_space spv (args(0).map_value (), msh);
  geopdes_space spp (args(1).map_value (), msh);
  Matrix coeff = args(3).matrix_value();

  SparseMatrix mat;

  if (!error_state)
    {

      const octave_idx_type nel = msh.nel (), nqn = msh.nqn (), ncomp = spv.ncomp (), ndof_spp = spp.ndof (), 
        nsh_max_spp = spp.nsh_max (), ndof_spv = spv.ndof (), nsh_max_spv = spv.nsh_max ();

#if OCTAVE_API_VERSION_NUMBER>37
      dim_vector dims (nel * nsh_max_spv * nsh_max_spp, 1);
      Array <octave_idx_type> I (dims, 0);
      octave_idx_type* Iptr = I.fortran_vec ();

      Array <octave_idx_type> J (dims, 0);
      octave_idx_type* Jptr = J.fortran_vec ();

      Array <double> V (dims, 0.0);    
      double* Vptr = V.fortran_vec ();
#else
      ColumnVector I (nel * nsh_max_spv * nsh_max_spp, 0);
      double* Iptr = I.fortran_vec ();

      ColumnVector J (nel * nsh_max_spv * nsh_max_spp, 0);
      double* Jptr = J.fortran_vec ();

      ColumnVector V (nel * nsh_max_spv * nsh_max_spp, 0.0);
      double* Vptr = V.fortran_vec ();
#endif

      octave_idx_type counter = 0, iel, inode, idof, jdof, icmp;
 
      {      
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

              double shpv[nsh_v][nqn][ncomp];
              double shgp[nsh_p][nqn][ncomp];
              octave_idx_type conn_v[nsh_max_spv];
              octave_idx_type conn_p[nsh_max_spp];

              for ( idof = 0; idof < nsh_p; idof++) 
                for ( inode = 0; inode < nqn; inode++)
                  for ( icmp = 0; icmp < ncomp; icmp++)
                    shgp[idof][inode][icmp] = spp.shape_function_gradients (0, icmp, inode, idof, iel);

              for ( jdof = 0; jdof < nsh_v; jdof++) 
                for ( inode = 0; inode < nqn; inode++)
                  for ( icmp = 0; icmp < ncomp; icmp++)
                    shpv[jdof][inode][icmp] = spv.shape_functions (icmp, inode, jdof, iel);

              // Cache element data to speed-up memory access
              spp.cache_element_connectivity (iel, (octave_idx_type*)conn_p);
              spv.cache_element_connectivity (iel, (octave_idx_type*)conn_v);


              for ( idof = 0; idof < nsh_p; idof++) 
                {
                  for ( jdof = 0; jdof < nsh_v; jdof++)
                    {
                      counter = jdof + nsh_v * (idof + nsh_p * iel);

                      Iptr[counter] = conn_p[idof] - 1;
                      Jptr[counter] = conn_v[jdof] - 1;
                      Vptr[counter] = 0.0;
                      for ( inode = 0; inode < nqn; inode++)
                        {
                          if (msh.weights (inode, iel) > 0.0)
                            {			
                              double s = 0.0;
                              for ( icmp = 0; icmp < ncomp; icmp++)
                                s += shgp[idof][inode][icmp] * shpv[jdof][inode][icmp];
                              Vptr[counter] += jacdet_weights[inode] * s;
                            }  
                        } // end for inode		  
                      counter++;
                    } // end for jdof
                } // end for idof
            } else {
            {warning_with_id ("geopdes:zero_measure_element", "op_v_gradp: element %d has 0 area (or volume)", iel);}
          }  // end for iel, if area > 0
      } // end of parallel section

      if (nargout == 1) 
        {
          mat = SparseMatrix (V, I, J, ndof_spp, ndof_spv, true);
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
