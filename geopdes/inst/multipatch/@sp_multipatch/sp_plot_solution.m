% SP_PLOT_SOLUTION: Plot the computed solution, given the degrees of freedom.
%
%   [eu, F] = sp_plot_solution (u, space, geometry, pts, [ncuts=2]);
%   [eu, F] = sp_plot_solution (u, space, geometry, [npts], [ncuts=2]);
%
% INPUT:
%     
%     u:           vector of dof weights
%     space:       object defining the discrete space (see sp_multipatch)
%     geometry:    an array of geometry structures (see mp_geo_load)
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

if (nargin < 4)
  npts = [];
end
if (nargin < 5)
  ncuts = [];
end

if (isa (space.sp_patch{1}, 'sp_vector'))
  disp ('Warning: a different scaling is used for each patch')
end

hold_flag = ishold ;
for iptc = 1:space.npatch
  if (isempty (space.dofs_ornt))
    sp_plot_solution (u(space.gnum{iptc}), space.sp_patch{iptc}, geometry(iptc), npts, ncuts);
  else
    sp_plot_solution (u(space.gnum{iptc}) .* space.dofs_ornt{iptc}.', space.sp_patch{iptc}, geometry(iptc), npts, ncuts);
  end
  hold on
end

if (~hold_flag)
  hold off
end

end