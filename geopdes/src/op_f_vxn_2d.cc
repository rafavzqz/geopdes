/* Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
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

DEFUN_DLD(op_f_vxn_2d, args, nargout,"\n\
OP_F_VXN_2D: assemble the vector r = [r(i)], with  r(i) = (f, v_i x n), with n the normal exterior vector. \n\
\n\
  mat = op_f_vxn_2d (spv, msh, coeff);\n\
\n\
INPUT:\n\
\n\
  spv:   structure representing the function space (see sp_vector/sp_eval_boundary_side)\n\
  msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_eval_boundary_side)\n\
  coeff: source function evaluated at the quadrature points\n\
\n\
OUTPUT:\n\
\n\
  mat: assembled right-hand side \n\
")
{
  
  octave_value_list retval;

  geopdes_mesh_normal msh (args(1).scalar_map_value ());
  geopdes_space       sp  (args(0).scalar_map_value (), msh);
  NDArray             coeff = args(2).array_value ();

  if (!error_state)
    {
      const octave_idx_type nel = msh.nel (), nqn = msh.nqn ();

     ColumnVector mat (sp.ndof (), 0.0);

      octave_idx_type iel, inode, idof;
      {
        double local_contribution;
        for (iel=0; iel < nel; iel++) 
          if (msh.area (iel) > 0)
            {
              const octave_idx_type nsh = sp.nsh (iel);
              double jacdet_weights[nqn];

            for (inode = 0; inode < nqn; inode++)
              jacdet_weights[inode] = msh.jacdet (inode, iel) *
                msh.weights (inode, iel) * coeff (inode, iel);
            
            double shp_x_n[nsh][nqn];
            octave_idx_type conn[nsh];

            for (idof = 0; idof < nsh; idof++) 
              for (inode = 0; inode < nqn; inode++)
                shp_x_n[idof][inode] = 
                  sp.shape_functions (0, inode, idof, iel)*
                  msh.normal (1, inode, iel) -
                  sp.shape_functions (1, inode, idof, iel)*
                  msh.normal (0, inode, iel);

            sp.cache_element_connectivity (iel, (octave_idx_type*)conn);

            for (idof = 0; idof < nsh; idof++) 
              for (inode = 0; inode < nqn; inode++)
                if (msh.weights (inode, iel) > 0.0)
                  {
                    local_contribution = jacdet_weights[inode] * 
                      shp_x_n[idof][inode];
                    mat(conn[idof]-1) += local_contribution;
                  }  // end for idof, for inode, if  
            } 
          else
            {
              {warning_with_id ("geopdes:zero_measure_element", "op_f_vxn_2d: element %d has 0 area", iel);}
            } // end for iel, if area > 0
      } // end of parallel region
      retval(0) = octave_value (mat);
    } // end if !error_state
  return retval;
}
