% MP_SP_DRCHLT_L2_PROJ: assign the degrees of freedom of Dirichlet boundaries through an L2 projection, in multipatch geometries.
%  This function from version 2 of GeoPDEs is deprecated.
%
%   [u, dofs] = mp_sp_drchlt_l2_proj (sp, msh, h, gnum, boundaries, refs)
%
% INPUT:
%
%  sp:         object representing the space of trial functions (see sp_bspline)
%  msh:        object containing the domain partition and the quadrature rule (see msh_cartesian)
%  h:          function handle to compute the Dirichlet condition
%  gnum:       global numbering of the degrees of freedom (see mp_interface)
%  boundaries: array of structures containing the information for the boundaries (see mp_geo_load)
%  refs:       boundary references on which a Dirichlet condition is imposed
%
% OUTPUT:
%
%  u:    assigned value to the degrees of freedom
%  dofs: global numbering of the corresponding basis functions
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
% Copyright (C) 2011, 2015 Rafael Vazquez
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

function [u, dofs] = mp_sp_drchlt_l2_proj (sp, msh, h, gnum, boundaries, refs)

  if (isa (sp, 'sp_multipatch'))
    warning ('For spaces of the class SP_MULTIPATCH, using the function SP_DRCHLT_L2_PROJ inside the class')
    [u, dofs] = sp_drchlt_l2_proj (sp, msh, h, refs);
    return
  end

  dofs = [];
  ndof = max ([gnum{:}]);
  M    = spalloc (ndof, ndof, ndof);
  rhs  = spalloc (ndof, 1, ndof);

  for iref = refs
    for bnd_side = 1:boundaries(iref).nsides
      iptc = boundaries(iref).patches(bnd_side);
      iside = boundaries(iref).faces(bnd_side);

      global_dofs = gnum{iptc}(sp{iptc}.boundary(iside).dofs);
      dofs = union (dofs, global_dofs);
% Restrict the function handle to the specified side, in any dimension, hside = @(x,y) h(x,y,iside)
      href = @(varargin) h(varargin{:},iref);
      f_one = @(varargin) ones (size(varargin{1}));

      M_side = op_u_v_tp (sp{iptc}.boundary(iside), sp{iptc}.boundary(iside), msh{iptc}.boundary(iside), f_one);
      M(global_dofs, global_dofs) = M(global_dofs, global_dofs) + M_side;

      rhs_side = op_f_v_tp (sp{iptc}.boundary(iside), msh{iptc}.boundary(iside), href);
      rhs(global_dofs) = rhs(global_dofs) + rhs_side;
    end
  end

  dofs = unique (dofs);
  u = M(dofs, dofs) \ full (rhs(dofs));

end
