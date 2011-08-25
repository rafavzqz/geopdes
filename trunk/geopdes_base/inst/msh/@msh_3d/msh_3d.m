% MSH_3D: constructor of the class for 3d tensor product meshes.
%
%     msh = msh_3d (breaks, qn, qw, geometry, 'option1', value1, ...);
%
% INPUTS:
%     
%     breaks:   breaks along each direction in parametric space (repetitions are ignored)
%     qn:       quadrature nodes along each direction in parametric space
%     qw:       quadrature weights along each direction in parametric space
%     geometry: structure representing the geometrical mapping
%    'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+-----------
%             boundary  |      true       |  compute the quadrature rule
%                       |                 |   also on the boundary
%   
% OUTPUT:
%
%     msh: object containing the following fields and methods
%
%     FIELD_NAME    (SIZE)                  DESCRIPTION
%     nel           (scalar)                number of elements of the partition
%     nel_dir       (1 x 3 vector)          number of elements in each parametric direction
%     nelcol        (scalar)                number of elements in one "column" of the mesh (actually, nel_dir(2)*nel_dir(3))
%     nqn           (scalar)                number of quadrature nodes per element
%     nqn_dir       (1 x 3 vector)          number of quadrature nodes per element in each parametric direction
%     breaks        (1 x 3 cell-array)      unique(breaks)
%     qn            (1 x 3 cell-array)      quadrature nodes along each direction in parametric domain
%     qw            (1 x 3 cell-array)      quadrature weights along each direction in parametric space
%     boundary      (1 x 6 struct-array)    it contains a one-dimensional 'msh' structure for each edge of the boundary (only when boundary is set to true)
%     map           (function handle)       a copy of the map handle of the geometry structure
%     map_der       (function handle)       a copy of the map_der handle of the geometry structure
%
%     METHOD NAME
%     msh_evaluate_col: computes the parameterization (and its derivatives) at
%                       the quadrature points in one column of the mesh, i.e.,
%                       fixing the element in the first parametric direction.
%     msh_eval_boundary_side: computes the parameterization in one boundary side
%                       of the domain.
%     msh_precompute:   computes, for all the mesh, any of the fields related
%                       to the quadrature rule, or all of them (except boundary)
%                       as in the mesh structure from previous versions.
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

function msh = msh_3d (breaks, qn, qw, geo, varargin)

  boundary = true;
  if (~isempty (varargin))
    if (~rem (length (varargin), 2) == 0)
      error ('msh_3d: options must be passed in the [option, value] format');
    end
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin {ii}, 'boundary'))
        boundary = varargin {ii+1};
      else
        error ('msh_3d: unknown option %s', varargin {ii});
      end
    end
  end

  msh.qn = qn;
  msh.qw = qw;

  % this only for precaution, breaks should already not have any repetitions
  msh.breaks = {unique(breaks{1}), unique(breaks{2}), unique(breaks{3})};

  msh.nel_dir(1) = numel (msh.breaks{1}) - 1;
  msh.nel_dir(2) = numel (msh.breaks{2}) - 1;
  msh.nel_dir(3) = numel (msh.breaks{3}) - 1;
  msh.nel = prod (msh.nel_dir);
  msh.nelcol = msh.nel_dir(2) * msh.nel_dir(3);

  qnu = qn{1};  qnv = qn{2}; qnw = qn{3};
  msh.nqn_dir(1) = size (qnu,1); 
  msh.nqn_dir(2) = size (qnv,1);
  msh.nqn_dir(3) = size (qnw,1);
  msh.nqn  = prod (msh.nqn_dir);
  
  if (boundary)
    for iside = 1:6
%%    ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2]
      ind = setdiff (1:3, ceil(iside/2)); 

      bndry.side_number = iside;
      bndry.qn = {qn{ind}};
      bndry.qw = {qw{ind}};

      bndry.nel_dir = msh.nel_dir(ind);
      bndry.nel = prod (bndry.nel_dir);
      bndry.nqn_dir = msh.nqn_dir(ind);
      bndry.nqn = prod (bndry.nqn_dir);

      msh.boundary(iside) = bndry;
    end
  else
    msh.boundary = [];    
  end

  msh.map = geo.map;
  msh.map_der = geo.map_der;

  msh.quad_nodes   = [];
  msh.quad_weights = [];
  msh.geo_map      = [];
  msh.geo_map_jac  = [];
  msh.jacdet       = [];

  msh = class (msh, 'msh_3d');

end
