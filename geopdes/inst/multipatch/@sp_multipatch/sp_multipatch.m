% SP_MULTIPATCH: Constructor of the class for multipatch spaces.
%
%     sp = sp_multipatch (spaces, msh, interfaces)
%     sp = sp_multipatch (spaces, msh, interfaces, boundary_interfaces)
%
% INPUTS:
%
%    spaces:     cell-array of space objects, one for each patch (see sp_scalar, sp_vector)
%    msh:        mesh object that defines the multipatch mesh (see msh_multipatch)
%    interfaces: information of connectivity between patches (see mp_geo_load)
%    boundary_interfaces: information of connectivity between boundary patches (see mp_geo_load)
%
% OUTPUT:
%
%    sp: object representing the discrete function space of vector-valued functions, with the following fields and methods:
%
%        FIELD_NAME      (SIZE)                       DESCRIPTION
%        npatch          (scalar)                      number of patches
%        ncomp           (scalar)                      number of components of the functions of the space (equal to msh.rdim)
%        ndof            (scalar)                      total number of degrees of freedom after gluing patches together
%        ndof_per_patch  (1 x npatch array)            number of degrees of freedom per patch
%        sp_patch        (1 x npatch cell-array)       the input spaces, one space object for each patch (see sp_scalar and sp_vector)
%        transform       (string)                      one of 'grad-preserving', 'curl-preserving' and 'div-preserving'
%        gnum            (1 x npatch cell-array)       global numbering of the degress of freedom (see mp_interface)
%        dofs_ornt       (1 x npatch cell-array)       global orientation of the degrees for freedom, for curl-conforming and div-conforming spaces
%        boundary        (1 x 1 object)                a (ndim-1) dimensional "sp_multipatch" object for the whole boundary 
%        dofs            (1 x ndof vector)             only for boundary spaces, degrees of freedom that do not vanish on the boundary
%        orientation     (1 x ndof vector)             only for boundary spaces with 'curl-preserving' transform, global orientation
%                                                       of the boundary dofs with respect to the volumetric one.
%        constructor     function handle               function handle to construct the same discrete space in a different msh
%
%       METHODS
%       Methods that give a structure with all the functions computed in a certain subset of the mesh
%         sp_evaluate_element_list: compute basis functions (and derivatives) in a given list of elements
%
%       Methods for post-processing, that require a computed vector of degrees of freedom
%         sp_h1_error:    compute the error in H1 norm
%         sp_l2_error:    compute the error in L2 norm
%         sp_hcurl_error: compute the error in H(curl) norm
%         sp_to_vtk:      export the computed solution to a pvd file, using a Cartesian grid of points on each patch
%
%       Methods for basic connectivity operations
%         sp_get_basis_functions: compute the functions that do not vanish in a given list of elements
%         sp_get_cells:           compute the cells on which a list of functions do not vanish
%         sp_get_neighbors:       compute the neighbors, functions that share at least one element with a given one
%
% Copyright (C) 2015 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function sp = sp_multipatch (spaces, msh, interfaces, boundary_interfaces)

  sp_class = class (spaces{1});
  if (~all (cellfun (@(x) isa (x, sp_class), spaces)))
    error ('All the spaces in the array should be of the same class')
  end
  aux = struct ([spaces{:}]);

  sp.npatch = numel (spaces);
  if (sp.npatch ~= msh.npatch)
    error ('The list of spaces does not correspond to the mesh')
  end


  sp.ncomp = spaces{1}.ncomp;
  sp.transform = spaces{1}.transform;
  
  if (~all ([aux.ncomp] == sp.ncomp))
    error ('The number of components should be the same for all the spaces')  
  end
  for iptc = 1:sp.npatch
    if (~strcmpi (spaces{iptc}.transform, sp.transform))
      error ('The transform to the physical domain should be the same for all the spaces')
    end
  end
  
  sp.ndof = 0;
  sp.ndof_per_patch = [aux.ndof];
  sp.sp_patch = spaces;

  if (strcmpi (sp_class, 'sp_scalar'))
    if (strcmpi (sp.transform, 'grad-preserving'))
      [sp.gnum, sp.ndof] = mp_interface (interfaces, spaces);
    elseif (strcmpi (sp.transform, 'integral-preserving'))
      sp.ndof = sum (sp.ndof_per_patch);
      Nf = cumsum ([0 sp.ndof_per_patch]);
      for iptc = 1:sp.npatch
        sp.gnum{iptc} = Nf(iptc)+1:Nf(iptc+1);
      end
    else
      error ('sp_multipatch: Unknown transformation')
    end
    sp.dofs_ornt = [];

  elseif (strcmpi (sp_class, 'sp_vector'))
    if (strcmpi (sp.transform, 'grad-preserving'))
      [sp.gnum, sp.ndof] = mp_interface_vector (interfaces, spaces);
      sp.dofs_ornt = [];
    elseif (strcmpi (sp.transform, 'curl-preserving'))
      [sp.gnum, sp.ndof, sp.dofs_ornt] = mp_interface_hcurl (interfaces, spaces);
    elseif (strcmpi (sp.transform, 'div-preserving'))
      [sp.gnum, sp.ndof, sp.dofs_ornt] = mp_interface_hdiv (interfaces, spaces, msh);
    else
      error ('sp_multipatch: Unknown transformation')
    end
  else
    error ('sp_multipatch: Unknown space class')
  end

  sp.interfaces = interfaces;
  
% Boundary construction
  if (nargin == 4 && ~isempty (msh.boundary) && ~isempty (spaces{1}.boundary))
    sp_bnd = cell (msh.boundary.npatch, 1);
    for iptc = 1:msh.boundary.npatch
      patch_number = msh.boundary.patch_numbers(iptc);
      side_number  = msh.boundary.side_numbers(iptc);
      sp_bnd{iptc} = spaces{patch_number}.boundary(side_number);
    end
    sp.boundary = sp_multipatch (sp_bnd, msh.boundary, boundary_interfaces);
    
    dofs = zeros (sp.boundary.ndof, 1);
    boundary_orientation = []; %boundary_orient = zeros (sp.boundary.ndof, 1);
    for iptc = 1:msh.boundary.npatch
      patch_number = msh.boundary.patch_numbers(iptc);
      side_number  = msh.boundary.side_numbers(iptc);
      dofs(sp.boundary.gnum{iptc}) = sp.gnum{patch_number}(sp.sp_patch{patch_number}.boundary(side_number).dofs);
      if (~isempty (sp.boundary.dofs_ornt))
        boundary_orientation(sp.boundary.gnum{iptc}) = sp.boundary.dofs_ornt{iptc} .* ...
          sp.dofs_ornt{patch_number}(sp.sp_patch{patch_number}.boundary(side_number).dofs);
      end
    end
    sp.boundary.dofs = dofs;
    sp.boundary.boundary_orientation = boundary_orientation;
    
  else
    sp.boundary = [];
  end
  
  sp.dofs = [];
  sp.boundary_orientation = [];
  
  sp.constructor = @(MSH) sp_multipatch (patches_constructor(spaces, MSH), MSH, interfaces);
    function spaux = patches_constructor (spaces, MSH)
      for ipatch = 1:MSH.npatch
        spaux{ipatch} = spaces{ipatch}.constructor(MSH.msh_patch{ipatch});
      end
    end

  sp = class (sp, 'sp_multipatch');
  
end
