% SP_REFINE: construct a refined space from a given one. The function only
%                refines the space, the mesh must be refined separately.
%
%     [sp_fine, Proj, Proj0, Proj1] = sp_refine (space, msh, nsub, [degree], [regularity]);
%
% The same number of subdivisions, degree and regularity is applied to every patch
%
% INPUTS:
%     
%     space:      the coarse space, an object of the sp_multipatch class (see sp_multipatch_C1)
%     msh:        an object of the msh_multipatch class, usuallly the refined mesh (see msh_multipatch)
%     nsub:       number of uniform subdivisions to apply on each knot span, and for each direction
%     degree:     degree of the fine space, and for each direction
%     regularity: regularity for the new space, and for each direction
%   
% OUTPUT:
%
%     sp_fine: the refined space, an object of the class sp_multipatch (see sp_multipatch)
%     Proj:    the coefficients relating 1D splines of the coarse and the fine spaces for each patch. 
%                A cell-array of size 1 x npatch, each entry containing the
%                coefficients for the patch (either for scalar or vector-valued spaces).
%     Proj0:   similar to Proj, for degree p and regularity r+1.
%     Proj1:   similar to Proj, for degree p-1 and regularity r.
%
% Copyright (C) 2015-2023 Rafael Vazquez
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

function [sp_fine, Proj, Proj0, Proj1] = sp_refine (space, msh, nsub, degree, regularity)

  if (nargin < 4)
    degree = space.sp_patch{1}.degree;
  end
  if (nargin < 5)
    regularity = degree - 2;
  end

  sp_fine = cell (1, space.npatch);
  Proj  = cell (1, space.npatch);
  Proj0 = cell (1, space.npatch);
  Proj1 = cell (1, space.npatch);
  for iptc = 1:space.npatch
    if (nargout > 1)
      [sp_fine{iptc}, Proj{iptc}] = sp_refine (space.sp_patch{iptc}, msh.msh_patch{iptc}, nsub, degree, regularity);
% Refine also the spaces of degree p, regularity r+1; and degree p-1, regularity r
      knots0_fine = kntrefine (space.knots0_patches{iptc}, nsub-1, degree, regularity+1);
      knots1_fine = kntrefine (space.knots1_patches{iptc}, nsub-1, degree-1, regularity);
      for idim = 1:msh.ndim
        Proj0{iptc}{idim} = basiskntins (degree(idim), space.knots0_patches{iptc}{idim}, knots0_fine{idim});
        Proj1{iptc}{idim} = basiskntins (degree(idim)-1, space.knots1_patches{iptc}{idim}, knots1_fine{idim});
      end
    else
      sp_fine{iptc} = sp_refine (space.sp_patch{iptc}, msh.msh_patch{iptc}, nsub, degree, regularity);
    end
  end
  sp_fine = sp_multipatch_C1 (sp_fine, msh, space.geometry, space.interfaces, space.vertices);

end