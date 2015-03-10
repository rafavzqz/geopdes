% SP_EVAL_BOUNDARY_SIDE: Construct the space structure of one side of the boundary.
%
%     sp_side = sp_eval_boundary_side (sp, msh_side)
%
% INPUTS:
%
%     sp:       space object (see sp_bspline_2d)
%     msh_side: mesh structure containing the information of the quadrature
%               rule on the boundary edge (see msh_2d/msh_eval_boundary_side)
%
% OUTPUT:
%
%     sp_side: structure that contains the following fields
%              (see the article for a detailed description)
%
%     FIELD_NAME      (SIZE)                                   DESCRIPTION
%     ncomp           (scalar)                                 number of components of the functions of the space (actually, 1)
%     nsh_max         (scalar)                                 maximum number of shape functions per element
%     nsh             (1 x msh_side.nel vector)                actual number of shape functions per each element
%     ndof            (scalar)                                 total number of degrees of freedom
%     connectivity    (nsh_max x msh_side.nel vector)          indices of basis functions that do not vanish in each element
%     shape_functions (msh_side.nqn x nsh_max x msh_side.nel)  basis functions evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010 Carlo de Falco
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

function sp_side = sp_eval_boundary_side (sp, msh_side)

  iside = msh_side.side_number;
  sp_side = sp.boundary(iside);

%%    ind  = [2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind  = [2 2 1 1] in 2D;
%%    ind2 = [1 1 2 2 3 3] in 3D,                  %ind2 = [1 1 2 2] in 2D
  ind2 = ceil (iside/2);
  ind = setdiff (1:msh_side.ndim+1, ind2);

  sp_univ = sp.sp_univ(ind);
  nsh_dim = {sp_univ.nsh};
  [nsh_grid{1:msh_side.ndim}] = ndgrid (nsh_dim{:});
  nsh = 1;
  for idim = 1:msh_side.ndim
    nsh = nsh .* nsh_grid{idim};
  end
  sp_side.nsh = nsh(:)';

  for idim = 1:msh_side.ndim
    csize = ones (1, 2*msh_side.ndim);
    csize([idim, msh_side.ndim+idim]) = [sp_univ(idim).nsh_max, msh_side.nel_dir(idim)];
    crep = [sp_univ.nsh_max, msh_side.nel_dir];
    crep([idim, msh_side.ndim+idim]) = 1;

    conn{idim} = reshape (sp_univ(idim).connectivity, csize);
    conn{idim} = repmat (conn{idim}, crep);
    conn{idim} = reshape (conn{idim}, [], msh_side.nel);
  end
  connectivity = zeros (sp_side.nsh_max, msh_side.nel);
  indices = ones (size (conn{1}));
  for idim = 1:msh_side.ndim
    indices = indices & conn{idim} ~= 0;
  end
  for idim = 1:msh_side.ndim
    conn{idim} = conn{idim}(indices);
  end
  connectivity(indices) = sub2ind ([sp_side.ndof_dir, 1], conn{:}); % The extra one makes things work in any dimension
  sp_side.connectivity = reshape (connectivity, sp_side.nsh_max, msh_side.nel);
  clear conn csize crep indices connectivity

  for idim = 1:msh_side.ndim
    ssize = ones (1, 3*msh_side.ndim);
    ssize([idim, msh_side.ndim+idim, 2*msh_side.ndim+idim]) = [msh_side.nqn_dir(idim), sp_univ(idim).nsh_max, msh_side.nel_dir(idim)];
    srep = [msh_side.nqn_dir, sp_univ.nsh_max, msh_side.nel_dir];
    srep([idim, msh_side.ndim+idim, 2*msh_side.ndim+idim]) = 1;
    shp{idim} = reshape (sp_univ(idim).shape_functions, ssize);
    shp{idim} = repmat (shp{idim}, srep);
    shp{idim} = reshape (shp{idim}, msh_side.nqn, sp_side.nsh_max, msh_side.nel);
  end
  sp_side.shape_functions = 1;
  for idim = 1:msh_side.ndim
    sp_side.shape_functions = sp_side.shape_functions .* shp{idim};
  end

% TODO: try to define the boundary as another spline space
%  and compute with sp_precompute
end
