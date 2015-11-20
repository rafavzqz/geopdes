% SP_REFINE: construct a refined space from a given one. The function only
%                refines the space, the mesh must be refined separately.
%
%     [sp_fine, Proj] = sp_refine (space, msh, nsub, degree, regularity);
%
% The same number of subdivision is applied to each component. The degree
%   and the regularity may change.
%
% INPUTS:
%     
%     space:      the coarse space, an object of the sp_scalar class (see sp_vector)
%     msh:        an object of the msh_cartesian class, usuallly the refined mesh (see msh_cartesian)
%     nsub:       number of uniform subdivisions to apply on each knot span, for each direction
%     degree:     degree of the fine space, for each component and for each direction
%     regularity: regularity for the new space, for each component and for each direction
%   
% OUTPUT:
%
%     sp_fine: the refined space, an object of the class sp_scalar (see sp_vector)
%     Proj:    the coefficients relating 1D splines of the coarse and the fine spaces
%               and for each component. Cell-array of size ncomp_(param) x ndim
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

function [sp_fine, Proj] = sp_refine (space, msh, nsub, degree, regularity)

% The number of scalar spaces can be equal to ncomp or to ncomp_param
  nspaces = numel (space.scalar_spaces);

  if (~iscell (degree))
    for icomp = 1:nspaces
      deg{icomp} = degree;
    end
  else
    deg = degree;
  end
  if (~iscell (regularity))
    for icomp = 1:nspaces
      reg{icomp} = regularity;
    end
  else
    reg = regularity;
  end

  scalar_spaces = cell (1, nspaces);
  Proj = cell (nspaces, msh.ndim);
  for icomp = 1:nspaces
    if (nargout == 2)
      [scalar_spaces{icomp}, Proj_comp] = sp_refine (space.scalar_spaces{icomp}, msh, nsub, deg{icomp}, reg{icomp});
      Proj(icomp,:) = Proj_comp(:);
    else
      scalar_spaces{icomp} = sp_refine (space.scalar_spaces{icomp}, msh, nsub, deg{icomp}, reg{icomp});
    end
  end
  sp_fine = sp_vector (scalar_spaces, msh, space.transform);
  
end