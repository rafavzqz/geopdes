/* Copyright (C) 2009 Carlo de Falco
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

DEFUN_DLD(op_gradu_gradv, args, nargout,"\n\
OP_GRADU_GRADV: assemble the stiffness matrix A = [a(i,j)], a(i,j) = (epsilon grad u_j, grad v_i).\n\
\n\
  mat = op_gradu_gradv (spu, spv, msh, epsilon);\n\
\n\
INPUT: \n\
\n\
  spu:     structure representing the space of trial functions (see sp_bspline_2d_phys) \n\
  spv:     structure representing the space of test functions  (see sp_bspline_2d_phys) \n\
  msh:     structure containing the domain partition and the quadrature rule (see msh_push_forward_2d) \n\
  epsilon: diffusion coefficient \n\
\n\
OUTPUT: \n\
\n\
  mat: assembled stiffness matrix \n\
")
{
  
  octave_value_list retval;

  geopdes_mesh  msh (args(2).map_value ());
  geopdes_space spu (args(0).map_value (), msh);
  geopdes_space spv (args(1).map_value (), msh);
  Matrix coeff = args(3).matrix_value();

  if (!error_state)
    {
      
      const octave_idx_type nel = msh.nel (), ndir = msh.ndir (), ncomp = spu.ncomp (), nqn = msh.nqn (), ndof_spu = spu.ndof (), nsh_max_spu = spu.nsh_max (), ndof_spv = spv.ndof (), nsh_max_spv = spv.nsh_max ();

      dim_vector dims (nel * nsh_max_spv * nsh_max_spu, 1);
      Array <octave_idx_type> I (dims, 0);
      Array <octave_idx_type> J (dims, 0);
      Array <double> V (dims, 0.0);
     
      SparseMatrix mat;

      octave_idx_type iel, inode, idof, jdof, icmp, idir;
        
      {
        octave_idx_type counter = 0;

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

              double shgv[nsh_v][nqn][ncomp][ndir];
              double shgu[nsh_u][nqn][ncomp][ndir];
              int conn_v[nsh_v];
              int conn_u[nsh_u];

              for ( idof = 0; idof < nsh_v; idof++) 
		{
                  for ( inode = 0; inode < nqn; inode++)
                    {
                      for ( icmp = 0; icmp < ncomp; icmp++)
                        {
                          for ( idir = 0; idir < ndir; idir++)
                            {
                              shgv[idof][inode][icmp][idir] = spv.shape_function_gradients (icmp, idir, inode, idof, iel);
                            }
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
                    for ( idir = 0; idir < ndir; idir++)
                      {
                        shgu[jdof][inode][icmp][idir] = spu.shape_function_gradients (icmp, idir, inode, jdof, iel);
                      }
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
                                for ( idir = 0; idir < ndir; idir++)
                                   s += shgv[idof][inode][icmp][idir] * 
                                     shgu[jdof][inode][icmp][idir];
  
                              V(counter) += s * jacdet_weights[inode];
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
            {warning_with_id ("geopdes:zero_measure_element", "op_gradu_gradv: element %d has 0 area (or volume)", iel);}
          }  // end for iel, if area > 0
      } // end of openmp parallel section

      if (nargout == 1) 
        {
          mat = SparseMatrix (V, I, J, ndof_spv, ndof_spu, true);
          retval(0) = octave_value (mat);
        } 
      else if (nargout == 3)
	{
          for ( icmp = 0; icmp < I.length(); icmp++) 
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

