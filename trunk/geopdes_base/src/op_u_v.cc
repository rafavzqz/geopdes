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

  const int nargin = args.length ();

  if (!error_state)
    {      

      octave_idx_type iel, giel, inode, idof, jdof, icmp;
      octave_idx_type nel = msh.nel ();
      Array <octave_idx_type> element_list;
      if (nargin > 4) 
        {
          element_list = args(4).octave_idx_type_vector_value ();
          nel = element_list.numel ();
        } 
      else 
        {
          dim_vector tmpd (nel, 1);
          element_list.resize (tmpd, 0);
          for (iel = 0; iel < nel; iel++)
            element_list (iel) = iel + 1;
        }

      const octave_idx_type ncomp = spu.ncomp (), nqn = msh.nqn (),
        ndof_spu = spu.ndof (), nsh_max_spu = spu.nsh_max (), ndof_spv = spv.ndof (), 
        nsh_max_spv = spv.nsh_max ();

      octave_idx_type nsh_u = spu.nsh_max ();
      octave_idx_type nsh_v = spv.nsh_max ();
      double jacdet_weights[nqn];
      double shpv[nsh_v][nqn][ncomp];
      double shpu[nsh_u][nqn][ncomp];
      int conn_v[nsh_v];
      int conn_u[nsh_u];

      dim_vector dims (nel * nsh_max_spv * nsh_max_spu, 1);
      Array <octave_idx_type> I (dims, 0);
      Array <octave_idx_type> J (dims, 0);
      Array <double> V (dims, 0.0);

      SparseMatrix mat;

#pragma omp parallel default (none) shared (msh, spu, spv, I, J, V, coeff)
      {
        octave_idx_type counter = 0;
#pragma omp for
        for (iel=0; iel < nel; iel++)
          {
            giel = element_list (iel) - 1;
            if (msh.area (iel) > 0.0)
              {

                for ( inode = 0; inode < nqn; inode++)
                  {
                    jacdet_weights[inode] = msh.jacdet (inode, giel) *
                      msh.weights (inode, giel) * coeff (inode, giel);
                  }

                // Cache element data to speed-up memory access
                spv.cache_element_shape_functions (iel, (double*)shpv);
                spu.cache_element_shape_functions (iel, (double*)shpu);

                for ( idof = 0; idof < nsh_v; idof++) 
		  conn_v[idof] = spv.connectivity (idof, iel) - 1;

                for ( jdof = 0; jdof < nsh_u; jdof++) 
                  conn_u[jdof] = spu.connectivity (jdof, iel) - 1;

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
                            if (msh.weights (inode, giel) > 0.0)
                              {
                                double s = 0.0;
                                for ( icmp = 0; icmp < ncomp; icmp++)
                                  s += shpv[idof][inode][icmp] * 
                                    shpu[jdof][inode][icmp];
                              
                                V(counter) += s * jacdet_weights[inode];

                              }  
                          } // end for inode
                      } // end for jdof
                  } // end for idof
              } else {
#pragma omp critical
              {warning_with_id ("geopdes:zero_measure_element", "op_u_v: element %d has 0 area (or volume)", giel);}
            }  // end for iel
          }  // end if area > 0
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
          retval(2) = octave_value (V);
          retval(1) = octave_value (J);
          retval(0) = octave_value (I);
        }

    } // end if !error_state
  return retval;
}


