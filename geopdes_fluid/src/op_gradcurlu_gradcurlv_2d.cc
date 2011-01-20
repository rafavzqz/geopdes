/* Copyright (C) 2010 Carlo de Falco

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

class geopdes_space_stream: public geopdes_space
{
public:
  geopdes_space_stream (const Octave_map& refsp, const geopdes_mesh& msh): geopdes_space (refsp, msh) {};
  inline double shape_function_grad_curl (octave_idx_type i, octave_idx_type j, octave_idx_type k, 
                                              octave_idx_type m, octave_idx_type n) const
  {
    return (shape_function_hessians (0, i, j, k, m, n) * pow (-1.0, i));
  }
};


DEFUN_DLD(op_gradcurlu_gradcurlv_2d, args, nargout,"")
{
  
  octave_value_list retval;

  geopdes_mesh         msh (args(2).map_value ());
  geopdes_space_stream spu (args(0).map_value (), msh);
  geopdes_space_stream spv (args(1).map_value (), msh);
  Matrix               mu = args(3).matrix_value();

  if (!error_state)
    {
      dim_vector dims (msh.nel () * spv.nsh_max () * spu.nsh_max (), 1);
      Array <octave_idx_type> I (dims, 0);
      Array <octave_idx_type> J (dims, 0);
      Array <double> V (dims, 0.0);
      
      SparseMatrix mat;

#pragma omp parallel default (none) shared (msh, spu, spv, I, J, V, mu, lambda)
      {
      octave_idx_type counter;

#pragma omp for
      for (octave_idx_type iel=0; iel < msh.nel (); iel++) 
        if (msh.area (iel) > 0.0)
          {
            for ( octave_idx_type idof(0); idof < spv.nsh (iel); idof++) 
              {
                for ( octave_idx_type jdof(0); jdof < spu.nsh (iel); jdof++) 
                  {

                    counter = jdof + spu.nsh (iel) * (idof + spv.nsh (iel) * iel);

                    I(counter) = spv.connectivity (idof, iel) - 1;
                    J(counter) = spu.connectivity (jdof, iel) - 1;
                    V(counter) = 0.0;
                    for ( octave_idx_type inode(0); inode < msh.nqn (); inode++) 
                      {
                        if (msh.weights (inode, iel) > 0.0)
                          {
                            double s = 0.0;
                            for (octave_idx_type icmp(0); icmp < 2; icmp++)
                              for (octave_idx_type jcmp(0); jcmp < 2; jcmp++)
                                s += mu(inode, iel) * 
                                  spu.shape_function_grad_curl (icmp, jcmp, inode, jdof, iel) *
                                  spv.shape_function_grad_curl (icmp, jcmp, inode, idof, iel);
                                                        
                            V(counter) += s * msh.jacdet (inode, iel) * msh.weights (inode, iel);	
                          }  
                      } // end for inode
                    //		  if (idof != jdof) // copy upper triangular part to lower
                    //		    { 
                    //		      I(counter) = J(counter-1);
                    //		      J(counter) = I(counter-1);
                    //		      V(counter) = V(counter-1);
                    //		      counter++;
                    //		    } 
                  } // end for jdof
              } // end for idof
          } else {
#pragma omp critical
          {warning_with_id ("geopdes:zero_measure_element", "op_su_ev: element %d has 0 measure", iel);}
        }// end for iel
      }// end of parallel region
      mat = SparseMatrix (V, I, J, spv.ndof (), spu.ndof (), true);
      retval(0) = octave_value (mat);
    } // end if !error_state
  return retval;
}

