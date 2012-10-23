/* Copyright (C) 2010, 2011 Carlo de Falco

   This program is free software; you can redistribute it and/or modify
   it under the terms of the GNU General Public License as published by
   the Free Software Foundation; either version 3 of the License, or
   (at your option) any later version.
 
   This program is distributed in the hope that it will be useful,
   but WITHOUT ANY WARRANTY; without even the implied warranty of
   MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
   GNU General Public License for more details.

   You should have received a copy of the GNU General Public License
   along with this program; if not, see <http://www.gnu.org/licenses/>.
*/


#include <iostream>
#include <octave/oct.h>
#include <octave/oct-map.h>
#include <lo-ieee.h>

// ABSTRACT BASE CLASSES

class geopdes_mesh_base
{
protected:
  octave_idx_type nqn_rep, nel_rep, ndir_rep;
public:
  octave_idx_type nqn  (void) const { return nqn_rep; }
  octave_idx_type nel  (void) const { return nel_rep; }
  octave_idx_type ndir (void) const { return ndir_rep; }
  virtual double jacdet  (octave_idx_type inode, octave_idx_type iel) const = 0;
  virtual double weights (octave_idx_type inode, octave_idx_type iel) const = 0;
  inline double area (octave_idx_type iel) const {
    double a = 0.0;
    for (octave_idx_type iqn(0); iqn < nqn (); iqn++)
      a += std::abs ((jacdet (iqn, iel)) * (weights (iqn, iel)));
    return a;
  }
  inline double volume (octave_idx_type iel) const { return area (iel); }
};

class geopdes_space_base
{
protected:
  octave_idx_type ndof_rep, nsh_max_rep, ncomp_rep;

public:
  octave_idx_type ndof    (void) const { return ndof_rep;    }
  octave_idx_type nsh_max (void) const { return nsh_max_rep; }
  octave_idx_type ncomp   (void) const { return ncomp_rep;   }

  virtual octave_idx_type nsh (octave_idx_type iel) const = 0; 
  virtual octave_idx_type connectivity (octave_idx_type ish, octave_idx_type iel) const = 0;
  virtual void cache_element_connectivity (octave_idx_type iel, octave_idx_type *cache) const
  {
    for (octave_idx_type ii = 0; ii < nsh (iel); ii++)
      cache [ii] = connectivity (ii, iel);
  };

  virtual double shape_functions (octave_idx_type i, octave_idx_type j, octave_idx_type k, octave_idx_type m) const = 0;
  virtual void cache_element_shape_functions (octave_idx_type iel, double *cache) const {  };
  

  virtual double shape_function_gradients (octave_idx_type i, octave_idx_type j, octave_idx_type k, octave_idx_type m, octave_idx_type n) const {return octave_NaN;};
  virtual void cache_element_shape_function_gradients (octave_idx_type iel, double *cache) const {  };

  virtual double shape_function_curls (octave_idx_type i, octave_idx_type j, octave_idx_type k, octave_idx_type m) const {return octave_NaN;};
  virtual void cache_element_shape_function_curls (octave_idx_type iel, double *cache) const { };

  virtual double shape_function_divs  (octave_idx_type i, octave_idx_type j, octave_idx_type k, octave_idx_type m) const {return octave_NaN;};
  virtual void cache_element_shape_function_divs (octave_idx_type iel, double *cache) const { };

};

// PRELOADED MSH AND SPACE CLASSES

class geopdes_mesh: public geopdes_mesh_base
{
protected:
  Matrix jacdet_rep, weights_rep;
  Octave_map msh;

public:
  geopdes_mesh (const Octave_map& refmsh)   
  { 
    msh = refmsh;

    nqn_rep  = msh.contents  ("nqn")(0).int_value ();
    nel_rep  = msh.contents  ("nel")(0).int_value ();
    ndir_rep = msh.contents  ("geo_map_jac")(0).array_value ().dim2 (); 
    jacdet_rep  = msh.contents ("jacdet")(0).matrix_value (); 
    weights_rep = msh.contents ("quad_weights")(0).matrix_value (); 
  }

  double jacdet  (octave_idx_type inode, octave_idx_type iel) const { return jacdet_rep  (inode, iel); }
  double weights (octave_idx_type inode, octave_idx_type iel) const { return weights_rep (inode, iel); }

};

class geopdes_mesh_normal: public geopdes_mesh
{
protected: 
  NDArray normal_rep;

public:
  geopdes_mesh_normal (const Octave_map& msh): geopdes_mesh (msh)
  {
    normal_rep = msh.contents ("normal")(0).array_value ();
  }

  double normal (octave_idx_type i, octave_idx_type inode, octave_idx_type iel) const {return normal_rep (i, inode, iel); }
};


class geopdes_space: public geopdes_space_base, protected geopdes_mesh
{
protected:

  double * shape_functions_rep, 
    * shape_function_gradients_rep, 
    * shape_function_curls_rep, 
    * shape_function_divs_rep, 
    * shape_function_hessians_rep;
  octave_idx_type *nsh_rep, *connectivity_rep;

  Array<octave_idx_type> nsh_rep_v, connectivity_rep_v;
  Octave_map sp;

public:
  geopdes_space (const Octave_map& refsp, const geopdes_mesh& msh): geopdes_mesh (msh) 
  { 
    sp = refsp;

    ndof_rep    = sp.contents ("ndof")(0).int_value (); 
    nsh_max_rep = sp.contents ("nsh_max")(0).int_value (); 
    ncomp_rep   = sp.contents ("ncomp")(0).int_value (); 

    nsh_rep_v           = (sp.contents ("nsh")(0).octave_idx_type_vector_value ()); 
    nsh_rep             = (octave_idx_type *) nsh_rep_v.data (); 

    connectivity_rep_v  = (sp.contents ("connectivity")(0).octave_idx_type_vector_value (false, false, true)); 
    connectivity_rep    = (octave_idx_type *) connectivity_rep_v.data ();

    shape_functions_rep = NULL;
    if (sp.contains ("shape_functions"))
      shape_functions_rep = (double *) (sp.contents ("shape_functions")(0).array_value ()).data (); 

    shape_function_gradients_rep = NULL;
    if (sp.contains ("shape_function_gradients")) 
      shape_function_gradients_rep = (double *) (sp.contents ("shape_function_gradients")(0).array_value ()).data (); 

    shape_function_curls_rep = NULL;
    if (sp.contains ("shape_function_curls"))
      shape_function_curls_rep = (double *) (sp.contents ("shape_function_curls")(0).array_value ()).data (); 

    shape_function_divs_rep = NULL;
    if (sp.contains ("shape_function_divs"))
      shape_function_divs_rep = (double *) (sp.contents ("shape_function_divs")(0).array_value ()).data ();         

    shape_function_hessians_rep = NULL;
    if (sp.contains ("shape_function_hessians"))
      shape_function_hessians_rep = (double *) (sp.contents ("shape_function_hessians")(0).array_value ()).data ();         

  }

  inline octave_idx_type nsh (octave_idx_type iel) const { return nsh_rep[iel]; }
  inline octave_idx_type connectivity (octave_idx_type ish, octave_idx_type iel) const { return connectivity_rep[ish + nsh_max () * iel]; }

  void cache_element_connectivity (octave_idx_type iel, octave_idx_type *cache) const {
    octave_idx_type s  = nsh_max ();
    octave_idx_type * start = & (connectivity_rep[s * iel]);
    memcpy (cache, start, s * sizeof (octave_idx_type));
  }

  inline double shape_functions (octave_idx_type i, octave_idx_type j, octave_idx_type k, octave_idx_type m) const {
    return shape_functions_rep[i + ncomp () * (j + nqn () * (k + nsh_max () * m))];
  }

  void cache_element_shape_functions (octave_idx_type iel, double *cache) const {
    octave_idx_type s = ncomp () * nqn () * nsh_max ();
    double * start = & (shape_functions_rep[s * iel]);
    memcpy (cache, start, s * sizeof (double));
  }

  inline double shape_function_curls (octave_idx_type i, octave_idx_type j, octave_idx_type k, octave_idx_type m) const {
    double res = (shape_function_curls_rep != NULL) ? shape_function_curls_rep [i + ncomp () * (j + nqn () * (k + nsh_max () * m))] : octave_NaN;
    return res;
  }
  
  inline double shape_function_curls (octave_idx_type i, octave_idx_type j, octave_idx_type k) const {
    double res = (shape_function_curls_rep != NULL) ? shape_function_curls_rep [i + nqn () * (j + nsh_max () * k)] : octave_NaN;
    return res;
  }


  inline double shape_function_divs (octave_idx_type j, octave_idx_type k, octave_idx_type m) const {
    double res = (shape_function_divs_rep != NULL) ? shape_function_divs_rep   [(j + nqn () * (k + nsh_max () * m))] : octave_NaN;
    return res;
  }
  
  void cache_element_shape_function_divs (octave_idx_type iel, double *cache) const {
    octave_idx_type s = nqn () * nsh_max ();
    double * start = & (shape_function_divs_rep[s * iel]);
    memcpy (cache, start, s * sizeof (double));
  }

  inline double shape_function_gradients (octave_idx_type i, octave_idx_type j, octave_idx_type k, octave_idx_type m, octave_idx_type n) const {
    double res =  (shape_function_gradients_rep != NULL) ? shape_function_gradients_rep[i + ncomp () * (j + ndir () * (k + nqn () * (m + nsh_max () * n)))] : octave_NaN;
    return res;
  }

  void cache_element_shape_function_gradients (octave_idx_type iel, double *cache) const {
    octave_idx_type s = ncomp () * ndir () * nqn () * nsh_max ();
    double * start = & (shape_function_gradients_rep[s * iel]);
    memcpy (cache, start, s * sizeof (double));
  }

  inline double shape_function_hessians (octave_idx_type c, octave_idx_type i, octave_idx_type j, octave_idx_type k, octave_idx_type m, octave_idx_type n) const {
    double res =  (shape_function_hessians_rep != NULL) ? shape_function_hessians_rep[c + ncomp () * (i + ndir () * (j + ndir () * (k + nqn () * (m + nsh_max () * n))))] : octave_NaN;
    return res;
  }
};
