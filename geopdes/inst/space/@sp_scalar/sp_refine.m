% SP_REFINE: construct a refined space from a given one. The function only
%                refines the space, the mesh must be refined separately.
%
%     [sp_fine, Proj] = sp_refine (space, msh, nsub, degree, regularity);
%
% INPUTS:
%     
%     space:      the coarse space, an object of the sp_scalar class (see sp_scalar)
%     msh:        an object of the msh_cartesian class, usuallly the refined mesh (see msh_cartesian)
%     nsub:       number of uniform subdivisions to apply on each knot span, and for each direction
%     degree:     degree of the fine space, and for each direction
%     regularity: regularity for the new space, and for each direction
%   
% OUTPUT:
%
%     sp_fine: the refined space, an object of the class sp_scalar (see sp_scalar)
%     Proj:    the coefficients relating 1D splines of the coarse and the fine spaces
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

  if (any (degree < space.degree))
    error ('The given degree should be greater or equal than the degree of the coarse space')
  end
  
  if (any (regularity > degree - 1))
    error ('The regularity cannot be greater than degree-1')
  end

  Proj = cell (1, msh.ndim);
  if (strcmpi (space.space_type, 'spline'))
    knots = kntrefine (space.knots, nsub-1, degree, regularity);
    sp_fine = sp_bspline (knots, degree, msh);

    if (nargout == 2)
      if (all (sp_fine.degree == space.degree))
        for idim = 1:msh.ndim
          knt_coarse = space.knots{idim};
          knt_fine = sp_fine.knots{idim};
          Proj{idim} = basiskntins (degree(idim), knt_coarse, knt_fine);
        end
      else
        error ('For the Proj matrix, the case of degree elevation is not implemented yet')
      end
    end
    
  elseif (strcmpi (space.space_type, 'nurbs'))
% Overkilling to compute the new weights, but does not change that much
    coefs = zeros ([4, size(space.weights)]);
    coefs(4,:,:,:) = space.weights;
    if (numel (space.knots) ~= 1)
      aux_nurbs = nrbdegelev (nrbmak (coefs, space.knots), degree - space.degree);
    else
      aux_nurbs = nrbdegelev (nrbmak (coefs, space.knots{1}), degree - space.degree);
    end
    
    [knots,~,new_knots] = kntrefine (aux_nurbs.knots, nsub-1, aux_nurbs.order-1, regularity);
    aux_nurbs = nrbkntins (aux_nurbs, new_knots);
    weights = reshape (aux_nurbs.coefs(4,:,:,:), [aux_nurbs.number, 1]); % The extra 1 makes things work in any dimension

    sp_fine = sp_nurbs (knots, degree, weights, msh);

    if (nargout == 2)
      if (all (sp_fine.degree == space.degree))
        for idim = 1:msh.ndim
          knt_coarse = space.knots{idim};
          knt_fine = sp_fine.knots{idim};
          Proj{idim} = basiskntins (degree(idim), knt_coarse, knt_fine);
        end
      else
        error ('For the Proj matrix, the case of degree elevation is not implemented yet')
      end
    end
  end

end
