% MP_SP_DRCHLT_L2_PROJ: assign the degrees of freedom of Dirichlet boundaries through an L2 projection, in multipatch geometries.
%
%   [u, dofs] = mp_sp_drchlt_l2_proj (sp, msh, h, boundaries, refs)
%
% INPUT:
%
%  sp:         object representing the multipatch space of trial functions (see sp_multipatch)
%  msh:        object containing the domain partition and the quadrature rule (see msh_multipatch)
%  h:          function handle to compute the Dirichlet condition
%  boundaries: array of structures containing the information for the boundaries (see mp_geo_load)
%  refs:       boundary references on which a Dirichlet condition is imposed
%
% OUTPUT:
%
%  u:    assigned value to the degrees of freedom
%  dofs: global numbering of the corresponding basis functions
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

function [u, dofs] = mp_sp_drchlt_l2_proj (space, msh, h, boundaries, refs, dummy)

  if (nargin == 6)
    error (['The function MP_SP_DRCHLT_L2_PROJ has changed in version 3, to work with multipatch classes. ' ...
        'The old version can be called with MP_SP_DRCHLT_L2_PROJ_OLD'])
  end

  M = spalloc (space.boundary.ndof, space.boundary.ndof, 3*space.boundary.ndof);
  rhs = zeros (space.boundary.ndof, 1);
  
  Nbnd = cumsum ([0, boundaries.nsides]);
  bnd_dofs = [];
  for iref = refs
    iref_patch_list = Nbnd(iref)+1:Nbnd(iref+1);
    href = @(varargin) h(varargin{:}, iref);
    f_one = @(varargin) ones (size(varargin{1}));
    
    M = M + op_u_v_mp (space.boundary, space.boundary, msh.boundary, f_one, iref_patch_list);
    rhs = rhs + op_f_v_mp (space.boundary, msh.boundary, href, iref_patch_list);

    boundary_gnum = space.boundary.gnum;
    bnd_dofs = union (bnd_dofs, [boundary_gnum{iref_patch_list}]);
  end
  
  u = M(bnd_dofs,bnd_dofs) \ rhs(bnd_dofs);
  dofs = space.boundary.dofs(bnd_dofs);

end
