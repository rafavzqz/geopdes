% MSH_2D: constructor of the class for 2d tensor product meshes.
%
%     msh = msh_2d (breaks, qn, qw, geometry, 'option1', value1, ...);
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
%               der2    |      false      |  compute second order derivatives
%                       |                 |   of the geometry at quad nodes
%   
% OUTPUT:
%
%     msh: object containing the following fields and methods
%
%     FIELD_NAME    (SIZE)                  DESCRIPTION
%     nel           (scalar)                number of elements of the partition
%     nel_dir       (1 x 2 vector)          number of elements in each parametric direction
%     nelcol        (scalar)                number of elements in one "column" of the mesh (actually, nel_dir(2))
%     nqn           (scalar)                number of quadrature nodes per element
%     nqn_dir       (1 x 2 vector)          number of quadrature nodes per element in each parametric direction
%     breaks        (1 x 2 cell-array)      unique(breaks)
%     qn            (1 x 2 cell-array)      quadrature nodes along each direction in parametric domain
%     qw            (1 x 2 cell-array)      quadrature weights along each direction in parametric space
%     boundary      (1 x 4 struct-array)    it contains a one-dimensional 'msh' structure for each edge of the boundary (only when boundary is set to true)
%     der2          (scalar)                an option to say whether the second derivative must also be computed (by default is set to false)
%     map           (function handle)       a copy of the map handle of the geometry structure
%     map_der       (function handle)       a copy of the map_der handle of the geometry structure
%     map_der2      (function handle)       a copy of the map_der2 handle of the geometry structure
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

function msh = msh_2d (breaks, qn, qw, geo, varargin)

  boundary = true;
  msh.der2     = false;
  if (~isempty (varargin))
    if (~rem (length (varargin), 2) == 0)
      error ('msh_2d: options must be passed in the [option, value] format');
    end
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin {ii}, 'der2'))
        msh.der2 = varargin {ii+1};
      elseif (strcmpi (varargin {ii}, 'boundary'))
        boundary = varargin {ii+1};
      else
        error ('msh_2d: unknown option %s', varargin {ii});
      end
    end
  end

  msh.qn = qn;
  msh.qw = qw;

  % this only for precaution, breaks should already not have any repetitions
  msh.breaks = {unique(breaks{1}), unique(breaks{2})};

  msh.nel_dir(1) = numel (msh.breaks{1}) - 1;
  msh.nel_dir(2) = numel (msh.breaks{2}) - 1;
  msh.nel = msh.nel_dir(1) * msh.nel_dir(2);
  msh.nelcol = msh.nel_dir(2);

  qnu = qn{1};  qnv = qn{2};
  msh.nqn_dir(1) = size (qn{1},1); 
  msh.nqn_dir(2) = size (qn{2},1);
  msh.nqn  = prod (msh.nqn_dir);
  
  if (boundary)
    for iside = 1:4
      ind = mod (floor ((iside+1)/2), 2) + 1;  %ind = [2 2 1 1];

      msh.boundary(iside).side_number = iside;

      msh.boundary(iside).breaks = msh.breaks{ind};
      msh.boundary(iside).nel = size (qn{ind},2);
      msh.boundary(iside).nqn = size (qn{ind},1);
    end
  else
    msh.boundary = [];    
  end

  msh.map = geo.map;
  msh.map_der = geo.map_der;

  if (isfield (geo, 'map_der2'))
    msh.map_der2 = geo.map_der2;
  else
    msh.map_der2 = [];
    if (msh.der2)
      warning ('msh_2d: a function to compute second order derivatives has not been provided')
    end
  end

  msh.quad_nodes   = [];
  msh.quad_weights = [];
  msh.geo_map      = [];
  msh.geo_map_jac  = [];
  msh.geo_map_der2 = [];
  msh.jacdet       = [];

  msh = class (msh, 'msh_2d');

end
