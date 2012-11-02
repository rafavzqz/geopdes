% SP_EVAL_MSH: Evaluate a function, given by its degrees of freedom, at the points given by a msh object.
%
%   [eu, F] = sp_eval_msh (u, space, msh);
%
% INPUT:
%     
%     u:         vector of dof weights
%     space:     object defining the discrete space (see sp_bspline_3d)
%     msh:       object defining the points where to evaluate (see msh_3d)
%
% OUTPUT:
%
%     eu: the function evaluated in the given points 
%     F:  grid points in the physical domain, that is, the mapped points
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
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

function [eu, F] = sp_eval_msh (u, space, msh)

  ndim = numel (msh.qn);

  F  = zeros (ndim, msh.nqn, msh.nel);
  eu = zeros (space.ncomp, msh.nqn, msh.nel);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col);

    uc_iel = zeros (size (sp_col.connectivity));
    uc_iel(sp_col.connectivity~=0) = ...
          u(sp_col.connectivity(sp_col.connectivity~=0));
    weight = repmat (reshape (uc_iel, [1, 1, sp_col.nsh_max, msh_col.nel]), ...
                                  [sp_col.ncomp, msh_col.nqn, 1, 1]);

    sp_col.shape_functions = reshape (sp_col.shape_functions, sp_col.ncomp, ...
                                      msh_col.nqn, sp_col.nsh_max, msh_col.nel);

    F(:,:,msh_col.elem_list) = msh_col.geo_map;
    eu(:,:,msh_col.elem_list) = reshape (sum (weight .* sp_col.shape_functions, 3), ...
                             sp_col.ncomp, msh_col.nqn, msh_col.nel);
  end

  if (space.ncomp == 1)
    eu = reshape (eu, msh.nqn, msh.nel);
  end

end
