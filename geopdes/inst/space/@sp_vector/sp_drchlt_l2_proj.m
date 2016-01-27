% SP_DRCHLT_L2_PROJ: assign the degrees of freedom of Dirichlet boundaries through an L2 projection.
%
%   [u, dofs] = sp_drchlt_l2_proj (sp, msh, h, sides)
%
% INPUT:
%
%  sp:    object defining the space of discrete functions (see sp_vector)
%  msh:   object defining the domain partition and the quadrature rule (see msh_cartesian)
%  h:     function handle to compute the Dirichlet condition
%  sides: boundary sides on which a Dirichlet condition is imposed
%
% OUTPUT:
%
%  u:    assigned value to the degrees of freedom
%  dofs: global numbering of the corresponding basis functions
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2010, 2011, 2015 Rafael Vazquez
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

function [u, dofs] = sp_drchlt_l2_proj (sp, msh, h, sides)

  rhs  = zeros (sp.ndof, 1);
  
  % In the 1D case, with an open knot vector,  it is not necessary to compute a projection.
  % For now it only works for scalars
  if (msh.ndim == 1)
    dofs = []; u = zeros (numel(sides), 1);
    for ii = 1:numel(sides)
      iside = sides(ii);
      dofs = union (dofs, sp.boundary(iside).dofs);
      if (iside == 1)
        u(ii) = h(msh.breaks{1}(1), iside);
      else
        u(ii) = h(msh.breaks{1}(end), iside); 
      end
    end
    u = u(:);
    return
  end

  dofs = [];
  nent = 0;
  for iside = sides
    nent = nent + msh.boundary(iside).nel * sp.boundary(iside).nsh_max^2;
    dofs = union (dofs, sp.boundary(iside).dofs);
  end

  rows = zeros (nent, 1);
  cols = zeros (nent, 1);
  vals = zeros (nent, 1);
  
  ncounter = 0;
  for iside = sides
% Restrict the function handle to the specified side, in any dimension, hside = @(x,y) h(x,y,iside)
    hside = @(varargin) h(varargin{:},iside);
    [rs, cs, vs] = op_u_v_tp (sp.boundary(iside), sp.boundary(iside), msh.boundary(iside));
    
    bnd_dofs = sp.boundary(iside).dofs;
    
    rows(ncounter+(1:numel(rs))) = bnd_dofs(rs);
    cols(ncounter+(1:numel(rs))) = bnd_dofs(cs);
    vals(ncounter+(1:numel(rs))) = vs;
    ncounter = ncounter + numel (rs);

    rhs(bnd_dofs) = rhs(bnd_dofs) + op_f_v_tp (sp.boundary(iside),msh.boundary(iside), hside);
  end

  M = sparse (rows(1:ncounter), cols(1:ncounter), vals(1:ncounter));
  u = M(dofs, dofs) \ rhs(dofs, 1);

end
