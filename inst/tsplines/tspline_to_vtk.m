% TSPLINE_TO_VTK: Export to VTK format for plotting.
%
% USAGE:
%
%  tspline_to_vtk (u, tspline, npts, filename, fieldname, [option])
%
% INPUT:
%     
%     u:         vector of dof weights
%     tspline:   Bezier extraction read from a file of the T-splines plug-in for Rhino
%                (see read_bezier_extraction)
%     npts:      number of points on each Bezier element, and for each parametric direction
%     filename:  name of the output file.
%     fieldname: how to name the saved variable in the vtk file
%     option:    accepted options are 'value' (default), 'gradient'
%
% OUTPUT:
%
%    none
%
% Copyright (C) 2009, 2010 Carlo de Falco
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

function tspline_to_vtk (u, tspline, npts, filename, fieldname, varargin)

  if (nargin == 5)
    option = 'value';
  else
    option = varargin{1};
  end

% Correct the dimension of npts
  if (numel (npts) == 1)
    npts = [npts, npts];
  elseif (numel (npts) > 2)
    npts = npts(1:2);
  end
  
  
% For plotting, the elements are not connected. 
%  We need at least 2 points per element in each parametric direction.
  npts = arrayfun (@(x) max(2,x), npts);
  
  [msh, space] = tspline_mesh_space (tspline, 'plot', npts);
  [eu, F] = tspline_sp_eval_msh (u, space, msh, option);

  if (space.ncomp == 1 || ~strcmp(option, 'gradient'))
    msh_to_vtu (F, eu, npts, filename, fieldname);
  else
    error ('tspline_to_vtk: For vector fields, the gradient cannot be saved')
  end

end
