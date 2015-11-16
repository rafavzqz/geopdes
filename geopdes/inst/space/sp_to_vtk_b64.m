% SP_TO_VTK: Export to VTK format for plotting (store data in
%            binary base64 encoded format).
%
%  sp_to_vtk (u, space, geometry, npts, filename,
%             fieldname, [option], [precision])
%  sp_to_vtk (u, space, geometry, pts, filename,
%             fieldname, [option], [precision])
%
% INPUT:
%     
%     u:          vector of dof weights
%     space:      object representing the space of discrete functions (see sp_scalar)
%     geometry:   geometry structure (see geo_load)
%     npts:       number of points along each parametric direction where to evaluate
%     pts:        cell array with the coordinates along each parametric direction of the points where to evaluate
%     filename:   name of the output file. 
%     fieldname:  how to name the saved variable in the vtk file
%     option:     accepted options are 'value' (default), 'gradient',
%                 and for vectors also 'curl', 'divergence'
%     precision:  accepted options are 32 (default), or 64, other
%                 values are silently ignored
%
%
% OUTPUT:
%
%    a vtk structured mesh file named <filename> is produced, data are
%    base64-encoded and appended at the end of the file to save disk
%    space and write time.
% 
% Copyright (C) 2009, 2010, 2014 Carlo de Falco
% Copyright (C) 2011, 2012 Rafael Vazquez
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

function sp_to_vtk_b64 (u, space, geometry, npts, filename, fieldname, varargin)

  if (nargin >= 6)
    option = 'value';
  else
    option = varargin{1};
  end

  if (~(nargin >= 7) || ~(precision == 64))
    precision = 32;
  end

  [eu, F] = sp_eval (u, space, geometry, npts, option);

  if (space.ncomp == 1 || ~strcmp(option, 'gradient'))
    msh_to_vtk_b64 (F, eu, filename, fieldname, precision);
  else
    error ('sp_to_vtk: For vector fields, the gradient cannot be saved')
  end

end
