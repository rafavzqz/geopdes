% EX_ARTICLEV3_EXAMPLE51: simple example to generate multipatch objects of a two-patch geometry.
%
% The geometry is given by the two squares, [0,1] x [0,1] and [1,2] x [0,1],
%  defined with biquadratic splines.
%
% The code generates the multipatch msh and space, as in Example 5.1 of the paper:
%
% R. Vazquez, A new design for the implementation of isogeometric analysis 
%   in Octave and Matlab: GeoPDEs 3.0. Tech, Report, IMATI-CNR, 2016
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

nrb(1) = nrbdegelev (nrb4surf ([0 0], [1 0], [0 1], [1 1]), [1 1]);
nrb(2) = nrbtform (nrb(1), vectrans([1 0 0]));

[geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (nrb);
npatch = numel (geometry);
for iptc = 1:npatch
  knots = geometry(iptc).nurbs.knots;
  [qn, qw] = msh_set_quad_nodes (knots, msh_gauss_nodes (geometry(iptc).nurbs.order));
  
  local_meshes{iptc} = msh_cartesian (knots, qn, qw, geometry(iptc), 'der2', false);
  local_spaces{iptc}  = sp_nurbs (geometry(iptc).nurbs, local_meshes{iptc});
end

msh = msh_multipatch (local_meshes, boundaries);
space = sp_multipatch (local_spaces, msh, interfaces, boundary_interfaces);

disp ('Local to global numbering in space.gnum')
disp(['PATCH 1: ', num2str(space.gnum{1}.')])
disp(['PATCH 2: ', num2str(space.gnum{2}.')])
fprintf('\n')

disp ('Boundary information in msh.boundary')
disp (['patch_numbers: ', num2str(msh.boundary.patch_numbers)])
disp (['side_numbers:  ', num2str(msh.boundary.side_numbers)])
fprintf('\n')

disp (['Number of boundary patches, in space.boundary.npatch: ', num2str(space.boundary.npatch)]);
fprintf('\n')
disp ('Local to global numbering in space.boundary.gnum:')
for ii = 1:space.boundary.npatch
  disp(['PATCH ',num2str(ii),':  ', num2str(space.boundary.gnum{ii}.')])  
end
fprintf('\n')

disp ('Boundary to volumetric numbering, in space.boundary.dofs')
disp (num2str(space.boundary.dofs.'))

fprintf('\n')
disp ('For the solution of problems in multipatch geometries, see also mp_solve_laplace.m')

%!test
%! nrb(1) = nrbdegelev (nrb4surf ([0 0], [1 0], [0 1], [1 1]), [1 1]);
%! nrb(2) = nrbtform (nrb(1), vectrans([1 0 0]));
%! [geometry, boundaries, interfaces, ~, boundary_interfaces] = mp_geo_load (nrb);
%! npatch = numel (geometry);
%! for iptc = 1:npatch
%!   knots = geometry(iptc).nurbs.knots;
%!   [qn, qw] = msh_set_quad_nodes (knots, msh_gauss_nodes (geometry(iptc).nurbs.order));
%!   local_meshes{iptc} = msh_cartesian (knots, qn, qw, geometry(iptc), 'der2', false);
%!   local_spaces{iptc}  = sp_nurbs (geometry(iptc).nurbs, local_meshes{iptc});
%! end
%! msh = msh_multipatch (local_meshes, boundaries);
%! space = sp_multipatch (local_spaces, msh, interfaces, boundary_interfaces);
%! assert (space.gnum{1}(:), [1 2 13 3 4 14 5 6 15].');
%! assert (space.gnum{2}(:), [13 7 8 14 9 10 15 11 12].');
%! assert (msh.boundary.patch_numbers, [1 1 1 2 2 2]);
%! assert (msh.boundary.side_numbers, [1 3 4 2 3 4]);
%! assert (space.boundary.npatch, 6)
%! assert (space.boundary.gnum{1}(:), [7 1 8].');
%! assert (space.boundary.gnum{2}(:), [7 2 9].');
%! assert (space.boundary.gnum{3}(:), [8 3 10].');
%! assert (space.boundary.gnum{4}(:), [11 4 12].');
%! assert (space.boundary.gnum{5}(:), [9 5 11].');
%! assert (space.boundary.gnum{6}(:), [10 6 12].');
%! assert (space.boundary.dofs(:), [3 2 6 10 7 11 1 5 13 15 8 12].');
