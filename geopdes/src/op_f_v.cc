/* Copyright (C) 2009, 2020 Carlo de Falco
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

template<class ARRAY, class NUMBER>
void op_f_v_worker (const geopdes_mesh &msh, const geopdes_space &sp,
                    const ARRAY &coeff, NUMBER *matptr)
{

  const octave_idx_type
    nel = msh.nel (),
    ncomp = sp.ncomp (), 
    nqn = msh.nqn (),
    nsh = sp.nsh_max ();

  double shp[nsh][nqn][ncomp];
  octave_idx_type conn[nsh];
  double jacdet_weights[nqn];
  
  octave_idx_type iel, inode, idof, icmp;
      
  for (iel=0; iel < nel; iel++) 
    if (msh.area (iel) > 0.0)
      {              
        for (inode = 0; inode < nqn; inode++)
          jacdet_weights[inode] = msh.jacdet (inode, iel) *
            msh.weights (inode, iel);
            
        sp.cache_element_shape_functions (iel, (double*)shp);
        sp.cache_element_connectivity (iel, (octave_idx_type*)conn);
            
        for (idof = 0; idof < nsh; idof++) 
          for (inode = 0; inode < nqn; inode++)
            if (msh.weights (inode, iel) > 0.0)
              {

                NUMBER s = 0.0;

                for (icmp = 0; icmp < ncomp; icmp++)
                  s += shp[idof][inode][icmp] * coeff (icmp, inode, iel);

                matptr[(conn[idof])-1] += jacdet_weights[inode] * s;
		  
              } // end for idof, for inode, if weight > 0
      }  // end for iel, if area > 0 
    else 
      {
        warning_with_id ("geopdes:zero_measure_element", "op_f_v: element %lld has 0 area (or volume)", static_cast<long long int> (iel));
      }
}

DEFUN_DLD(op_f_v, args, nargout,"\n\
OP_F_V: assemble the right-hand side vector r = [r(i)], with  r(i) = (f, v_i). \n\
\n\
  mat = op_f_v (spv, msh, coeff); \n\
\n\
INPUT:\n\
\n\
  spv:   structure representing the function space (see sp_scalar/sp_evaluate_col) \n\
  msh:   structure containing the domain partition and the quadrature rule (see msh_cartesian/msh_evaluate_col) \n\
  coeff: source function evaluated at the quadrature points \n\
\n\
OUTPUT:\n\
\n\
  mat: assembled right-hand side \n\
")
{
  
  octave_value_list retval;

  if (nargout != 1 || args.length () != 3)
    print_usage ();
    
  geopdes_mesh  msh (args(1).scalar_map_value ());
  geopdes_space sp  (args(0).scalar_map_value (), msh);
  dim_vector    idx (sp.ncomp (), msh.nqn (), msh.nel ());
  

  if (args(2).iscomplex ())
    {
      Array<Complex> coeff = args(2).complex_array_value ().reshape (idx);
      ComplexColumnVector mat (sp.ndof (), 0.0);
      Complex *matptr = mat.fortran_vec ();
      op_f_v_worker (msh, sp, coeff, matptr);
      retval(0) = octave_value (mat);
    }
  else
    {
      Array<double> coeff = args(2).array_value ().reshape (idx);
      ColumnVector  mat (sp.ndof (), 0.0);
      double *matptr = mat.fortran_vec ();
      op_f_v_worker (msh, sp, coeff, matptr);
      retval(0) = octave_value (mat);
    }
 

  return retval;
}

