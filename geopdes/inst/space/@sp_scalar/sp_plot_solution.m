% SP_PLOT_SOLUTION: Plot the computed solution, given the degrees of freedom.
%
%   [eu, F] = sp_plot_solution (u, space, geometry, pts, [ncuts=2]);
%   [eu, F] = sp_plot_solution (u, space, geometry, [npts=51], [ncuts=2]);
%
% INPUT:
%     
%     u:           vector of dof weights
%     space:       object defining the discrete space (see sp_scalar)
%     geometry:    geometry structure (see geo_load)
%     pts:         cell array with coordinates of points along each parametric direction
%     npts:        number of points along each parametric direction
%     ncuts:       only for volumetric domains, number of internal "cuts" in each parametric direction.
%                    The zero value will plot the solution on the boundary.
%
%    This function only plots the value of the solution. To plot other
%     quantities, such as the gradient, compute them with sp_eval.
%
% Copyright (C) 2016 Rafael Vazquez
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

function sp_plot_solution (u, space, geometry, npts, ncuts)

ndim = numel (space.knots);

if (nargin < 4 || isempty (npts))
  npts = 51 * ones (1, ndim);
else
  if (iscell (npts))
    vtk_pts = npts;
    npts = cellfun (@numel, vtk_pts);
  elseif (numel (npts) == 1)
    npts = npts * ones (1, ndim);
  end
end

if (ndim == 3)
  if (nargin < 5 || isempty (ncuts))
    ncuts = 2 * ones (1, ndim);
  elseif (numel (ncuts) == 1)
    ncuts = ncuts * ones (1, ndim);
  end
end

degree = space.degree;
if (~exist ('vtk_pts', 'var'))
  for idim = 1:ndim
    vtk_pts{idim} = linspace (space.knots{idim}(1+degree(idim)), space.knots{idim}(end-degree(idim)), npts(idim));
  end
end

[eu, F] = sp_eval (u, space, geometry, vtk_pts);
rdim = size (F, 1);

if (ndim == 1)
  if (rdim == 1)
    plot (F, eu)
  elseif (rdim == 2)
    plot3 (F(1,:), F(2,:), eu)
  else
    error ('The plot on 3D curves has not been implemented yet')
  end
elseif (ndim == 2)
  if (rdim == 2)
    [X, Y]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)));
    surf (X, Y, eu)
  elseif (rdim == 3)
    [X, Y, Z]  = deal (squeeze(F(1,:,:)), squeeze(F(2,:,:)), squeeze(F(3,:,:)));
    surf (X, Y, Z, eu)
  end
elseif (ndim == 3)
  hold_flag = ishold;

  for idim = 1:ndim
    plot_pts = vtk_pts;
    plot_pts{idim} = linspace(space.knots{idim}(1+degree(idim)), space.knots{idim}(end-degree(idim)), ncuts(idim)+2);
    [eu, F] = sp_eval (u, space, geometry, plot_pts);
    indices = {1:npts(1), 1:npts(2), 1:npts(3)};
    
    if (ncuts(idim) > 0)
      cuts = 2:ncuts(idim)+1;
    else
      cuts = 1:ncuts(idim)+2;
    end
    for ii = cuts
      indices{idim} = ii;
      surf (squeeze(F(1,indices{:})), squeeze(F(2,indices{:})), squeeze(F(3,indices{:})), squeeze(eu(indices{:})));
      hold on
    end
  end

  if (~hold_flag)
    hold off;
  end
end

end