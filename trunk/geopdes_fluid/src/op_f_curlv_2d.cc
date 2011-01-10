/* Copyright (C) 2009 Carlo de Falco

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

DEFUN_DLD(op_f_curlv_2d, args, nargout,"")
{
  
  octave_value_list retval;

  geopdes_mesh  msh (args(1).map_value ());
  geopdes_space sp  (args(0).map_value (), msh);
  dim_vector    idx (2, msh.nqn (), msh.nel ());
  NDArray       coeff = args(2).array_value ().reshape (idx);

  if (!error_state)
    {

      ColumnVector mat (sp.ndof (), 0.0);

#pragma omp parallel default (none) shared (msh, sp, idx, mat, coeff)
      {
        double local_contribution;

#pragma omp for
        for (octave_idx_type iel=0; iel < msh.nel (); iel++) 
          if (msh.area (iel) > 0.0)
            {
              for (octave_idx_type idof(0); idof < sp.nsh (iel); idof++) 
                {
                  for ( octave_idx_type inode(0); inode < msh.nqn (); inode++) 
                    {
                      if (msh.weights (inode, iel) > 0.0)
                        {
                          double s = 0.0;
                          for (octave_idx_type icmp(0); icmp < 2; icmp++)
                            s += pow (-1.0, icmp) * sp.shape_function_gradients (0, 1-icmp, inode, idof, iel) * coeff (icmp, inode, iel);
                          local_contribution = msh.jacdet  (inode, iel) * msh.weights (inode, iel) * s;
#pragma omp critical
                          mat(sp.connectivity (idof, iel) - 1) += local_contribution;			  
                        }  
                    } // end for inode
                } // end for idof
            } else {
#pragma omp critical
            {warning_with_id ("geopdes:zero_measure_element", "op_f_v: element %d has 0 area (or volume)", iel);}
          }  // end for iel, if area > 0
      }  // end for parallel region
      retval(0) = octave_value (mat);
    } // end if !error_state
  return retval;
}

