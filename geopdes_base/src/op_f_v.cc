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

DEFUN_DLD(op_f_v, args, nargout,"\n\
OP_F_V: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i). \n\
\n\
  mat = op_f_v (spv, msh, coeff); \n\
\n\
INPUT:\n\
\n\
  spv:   structure representing the function space (see sp_bspline_2d_phys) \n\
  msh:   structure containing the domain partition and the quadrature rule (see msh_push_forward_2d) \n\
  coeff: source function evaluated at the quadrature points \n\
\n\
OUTPUT:\n\
\n\
  mat: assembled right-hand side \n\
")
{
  
  octave_value_list retval;

  geopdes_mesh  msh (args(1).map_value ());
  geopdes_space sp  (args(0).map_value (), msh);
  dim_vector idx (sp.ncomp (), msh.nqn (), msh.nel ());
  NDArray  coeff = args(2).array_value ().reshape (idx);

  if (!error_state)
    {
      ColumnVector mat (sp.ndof (), 0.0);
      for (octave_idx_type iel(0); iel < msh.nel (); iel++) 
        if (msh.area (iel) > 0.0)
	  {
	    for (octave_idx_type idof(0); idof < sp.nsh (iel); idof++) 
	      {
                for ( octave_idx_type inode(0); inode < msh.nqn (); inode++) 
                  {
                    if (msh.weights (inode, iel) > 0.0)
                      {
                        double s = 0.0;
                        for (octave_idx_type icmp(0); icmp < sp.ncomp (); icmp++)
                            s += sp.shape_functions (icmp, inode, idof, iel) * coeff (icmp, inode, iel);

                        mat(sp.connectivity (idof, iel) - 1) += 
                          msh.jacdet  (inode, iel) *
                          msh.weights (inode, iel) *
                          s;			  
                      }  
                  } // end for inode
	      } // end for idof
          } else {
          warning_with_id ("geopdes:zero_measure_element", "op_f_v: element %d has 0 area (or volume)", iel);
        }  // end for iel, if area > 0
      retval(0) = octave_value (mat);
    } // end if !error_state
  return retval;
}

