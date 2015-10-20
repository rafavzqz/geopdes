% MSH_BOUNDARY_SIDE_FROM_INTERIOR: create a mesh object with quadrature points on one boundary side, 
%    but taking into account the information from the interior. This is necessary when imposing 
%    Dirichlet conditions in a weak form, or to compute normal derivatives, for instance.
%
%     msh_side = msh_boundary_side_from_interior (msh, iside);
%
% INPUTS:
%     
%    msh:   mesh object (see msh_cartesian)
%    iside: number of the boundary side to compute, from 1 to 2*msh.ndim (see the file geo_specs for documentation about face numbering)
%
% OUTPUT:
%
%     msh_side: mesh object, with quadrature points only on the chosen side (see msh_cartesian)
%
% Copyright (C) 2014, 2015 Rafael Vazquez
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

function msh_side_from_interior = msh_boundary_side_from_interior (msh, iside)

  brk_bnd = msh.breaks; qn_bnd = msh.qn; qw_bnd = msh.qw;
  ind2 = ceil (iside/2);
  if (mod (iside, 2) == 1)
    brk_bnd{ind2} = brk_bnd{ind2}(1:2);
    qn_bnd{ind2} = brk_bnd{ind2}(1);
    qw_bnd{ind2} = 1;
  else
    brk_bnd{ind2} = brk_bnd{ind2}(end-1:end);
    qn_bnd{ind2} = brk_bnd{ind2}(end);
    qw_bnd{ind2} = 1;
  end
  
  geo.map = msh.map; geo.map_der = msh.map_der; geo.map_der2 = msh.map_der2;
  geo.rdim = msh.rdim;

  on_off = warning ('query', 'geopdes:check_quadrature');
  warning ('off', 'geopdes:check_quadrature');
  msh_side_from_interior = msh_cartesian (brk_bnd, qn_bnd, qw_bnd, geo, 'boundary', false);
  warning (on_off.state, 'geopdes:check_quadrature');
%   msh_side_from_interior = msh_precompute (msh_side_from_interior);

end
