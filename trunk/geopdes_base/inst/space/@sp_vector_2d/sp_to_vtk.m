% SP_TO_VTK: Export to VTK format for plotting.
%
%  sp_to_vtk (u, space, geometry, npts, filename, fieldname)
%  sp_to_vtk (u, space, geometry, {upts, vpts}, filename, fieldname)
%
% INPUT:
%     
%     u:          vector of dof weights
%     space:      structure representing the space of discrete functions (see sp_bspline_2d_phys)
%     geometry:   geometry structure (see geo_load)
%     npts:       number of points along each parametric direction where to evaluate
%     upts, vpts: points along each parametric direction where to evaluate
%     filename:   name of the output file. 
%     fieldname:  how to name the saved variable in the vtk file
%
% OUTPUT:
%
%    none    
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

function sp_to_vtk (u, space, geometry, npts, filename, fieldname)

  if (iscell (npts))
    pts = npts;
  elseif (isvector (npts))
    pts = {(linspace (0, 1, npts(1))), (linspace (0, 1, npts(2)))};
  end

  [p, F] = sp_eval (u, space, geometry, pts);

  msh_to_vtk (F, p, filename, fieldname);

end
