% SP_TO_VTK: Export to VTK format for plotting.
%
%  sp_to_vtk (u, space, geometry, npts, filename, fieldname, [options], [lambda_lame, mu_lame])
%  sp_to_vtk (u, space, geometry, pts, filename, fieldname, [options], [lambda_lame, mu_lame])
%
% INPUT:
%     
%     u:           vector of dof weights
%     space:       object representing the space of discrete functions (see sp_vector)
%     geometry:    geometry structure (see geo_load)
%     npts:        number of points along each parametric direction where to evaluate
%     pts:         cell array with the coordinates along each parametric direction of the points where to evaluate
%     filename:    name of the output file. 
%     fieldnames:  how to name the saved variables in the vtk file
%     options:     cell array with the fields to plot
%                   accepted options are 'value' (default), 'gradient',
%                   and for vectors also 'curl', 'divergence', 'stress'
%     lambda_lame: function handle to the first Lame coefficient (only needed to compute 'stress')
%     mu_lame:     function handle for the second Lame coefficient (only needed to compute 'stress')
%
% OUTPUT:
%
%    none    
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2012, 2015 Rafael Vazquez
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

function sp_to_vtk (u, space, geometry, npts, filename, fieldname, varargin)

  [eu, F] = sp_eval (u, space, geometry, npts, varargin{:});

  msh_to_vtk (F, eu, filename, fieldname);

end
