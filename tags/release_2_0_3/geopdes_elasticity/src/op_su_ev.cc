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

class geopdes_space_elasticity: public geopdes_space
{
public:
  geopdes_space_elasticity (const Octave_map& refsp, const geopdes_mesh& msh): geopdes_space (refsp, msh) {};
  inline double shape_function_strain_tensor (octave_idx_type i, octave_idx_type j, octave_idx_type k, 
                                              octave_idx_type m, octave_idx_type n) const
  {
    return ((shape_function_gradients (i, j, k, m, n) + shape_function_gradients (j, i, k, m, n))/2.0);
  }
};


DEFUN_DLD(op_su_ev, args, nargout,"\n\
 OP_SU_EV: assemble the matrix A = [a(i,j)], a(i,j) = 1/2 (sigma (u_j), epsilon (v_i)).\n\
\n\
   mat = op_su_ev (spu, spv, msh, lambda, mu);\n\
\n\
 INPUT:\n\
    \n\
   spu:        structure representing the space of trial functions (see sp_bspline_2d_phys)\n\
   spv:        structure representing the space of test functions  (see sp_bspline_2d_phys)\n\
   msh:        structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)\n\
   lambda, mu: Lame' coefficients\n\
\n\
 OUTPUT:\n\
\n\
   mat: assembled matrix\n\
 \n\
 Copyright (C) 2009, 2010, 2011 Carlo de Falco\n\
 Copyright (C) 2011 Rafael Vazquez\n\
\n\
    This program is free software: you can redistribute it and/or modify\n\
    it under the terms of the GNU General Public License as published by\n\
    the Free Software Foundation, either version 3 of the License, or\n\
    (at your option) any later version.\n\
")
{
  
  octave_value_list retval;

  geopdes_mesh             msh (args(2).map_value ());
  geopdes_space_elasticity spu (args(0).map_value (), msh);
  geopdes_space_elasticity spv (args(1).map_value (), msh);
  Matrix                   lambda = args(3).matrix_value();
  Matrix                   mu     = args(4).matrix_value();

  if (!error_state)
    {
      const octave_idx_type nel = msh.nel (), ncomp = spu.ncomp (), 
        nqn = msh.nqn (), ndof_spu = spu.ndof (), nsh_u = spu.nsh_max (), 
        ndof_spv = spv.ndof (), nsh_v = spv.nsh_max ();

      double jacdet_weights[nqn];
      double coeff_mu[nqn];
      double coeff_lambda[nqn];
      double shgv[nsh_v][nqn][ncomp][ncomp];
      double shgu[nsh_u][nqn][ncomp][ncomp];
      double shdivv[nsh_v][nqn];
      double shdivu[nsh_u][nqn];
      int conn_v[nsh_v];
      int conn_u[nsh_u];

      dim_vector dims (nel * nsh_v * nsh_u, 1);

#if OCTAVE_API_VERSION_NUMBER>37
      Array <octave_idx_type> I (dims, 0);
      octave_idx_type* Iptr = I.fortran_vec ();

      Array <octave_idx_type> J (dims, 0);
      octave_idx_type* Jptr = J.fortran_vec ();

      Array <double> V (dims, 0.0);
      double* Vptr = V.fortran_vec ();
#else
      ColumnVector I (nel * nsh_v * nsh_u, 0);
      double* Iptr = I.fortran_vec ();

      ColumnVector J (nel * nsh_v * nsh_u, 0);
      double* Jptr = J.fortran_vec ();

      ColumnVector V (nel * nsh_v * nsh_u, 0.0);
      double* Vptr = V.fortran_vec ();
#endif
      
      SparseMatrix mat;

      octave_idx_type counter = 0, iel, inode, idof, jdof, icmp, jcmp;
      for (iel=0; iel < nel; iel++) 
        if (msh.area (iel) > 0.0)
          {
             
            for (inode = 0; inode < nqn; inode++)
              {
                jacdet_weights[inode] = msh.jacdet (inode, iel) *
                  msh.weights (inode, iel);
                coeff_mu[inode] = mu (inode, iel);
                coeff_lambda[inode] = lambda (inode, iel);
              }

            // Cache element data to speed-up memory access
            spv.cache_element_shape_function_divs (iel, (double*)shdivv);
            spu.cache_element_shape_function_divs (iel, (double*)shdivu);

            spu.cache_element_connectivity (iel, (octave_idx_type*)conn_u);
            spv.cache_element_connectivity (iel, (octave_idx_type*)conn_v);


            for (idof = 0; idof < spv.nsh (iel); idof++) 
              for (inode = 0; inode < nqn; inode++)
                for (icmp = 0; icmp < ncomp; icmp++)
                  for (jcmp = 0; jcmp < ncomp; jcmp++)
                    shgv[idof][inode][icmp][jcmp] = spv.shape_function_strain_tensor (icmp, jcmp, inode, idof, iel);
            

            for (jdof = 0; jdof < spu.nsh (iel); jdof++) 
              for (inode = 0; inode < nqn; inode++)
                for (icmp = 0; icmp < ncomp; icmp++)
                  for (jcmp = 0; jcmp < ncomp; jcmp++)
                    shgu[jdof][inode][icmp][jcmp] = spu.shape_function_strain_tensor (icmp, jcmp, inode, jdof, iel);
            

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
                          for (jcmp = 0; jcmp < ncomp; jcmp++)
                            s += shgu[jdof][inode][icmp][jcmp] * shgv[idof][inode][icmp][jcmp];
                                                        
                        Vptr[counter] += jacdet_weights[inode] * 
                          (2 * s * coeff_mu[inode] + coeff_lambda[inode] * 
                           shdivu[jdof][inode] * shdivv[idof][inode]);
                      }   // end for inode, if
                } // end for idof, for jdof
          } 
        else 
          {
            warning_with_id ("geopdes:zero_measure_element", "op_su_ev: element %d has 0 measure", iel);
          }// end for iel

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

