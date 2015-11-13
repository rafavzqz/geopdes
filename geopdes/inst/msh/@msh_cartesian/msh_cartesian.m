% MSH_CARTESIAN: constructor of the msh class for Cartesian meshes.
%
%     msh = msh_cartesian (breaks, qn, qw, geometry, 'option1', value1, ...);
%
% INPUTS:
%     
%     breaks:   breaks along each direction in parametric space (repetitions are ignored)
%     qn:       quadrature nodes along each direction in parametric space. For direction j, the size must be nqn_dir(j) x nel_dir(j)
%     qw:       quadrature weights along each direction in parametric space
%     geometry: structure representing the geometrical mapping
%    'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value     |  Meaning
%           ------------+---------------------+-----------
%             boundary  |      true           | compute the quadrature rule
%                       |                     |   also on the boundary
%             der2      | depends on geometry | compute second derivatives
%   
% OUTPUT:
%
%     msh: object containing the following fields and methods
%
%     FIELD_NAME    (SIZE)                    DESCRIPTION
%     ndim          (scalar)                  dimension of the parametric space
%     rdim          (scalar)                  dimension of the physical space
%     nel           (scalar)                  number of elements of the partition
%     nel_dir       (1 x ndim vector)         number of elements in each parametric direction
%     nqn           (scalar)                  number of quadrature nodes per element
%     nqn_dir       (1 x ndim vector)         number of quadrature nodes per element in each parametric direction
%     breaks        (1 x ndim cell-array)     unique(breaks)
%     qn            (1 x ndim cell-array)     quadrature nodes along each direction in parametric domain
%     qw            (1 x ndim cell-array)     quadrature weights along each direction in parametric space
%     boundary      (1 x 2*ndim struct-array) it contains a (ndim-1)-dimensional 'msh_cartesian' object for each side of the boundary (only when boundary is set to true)
%     map           (function handle)         a copy of the map handle of the geometry structure
%     map_der       (function handle)         a copy of the map_der handle of the geometry structure
%     map_der2      (function handle)         a copy of the map_der2 handle of the geometry structure
%
%     METHOD NAME
%     msh_evaluate_col: compute the parameterization (and its derivatives) at
%                       the quadrature points in one column of the mesh, i.e.,
%                       fixing the element in the first parametric direction.
%     msh_eval_boundary_side: compute the parameterization in one boundary side
%                       of the domain.
%     msh_precompute:   compute, for the entire mesh, any of the fields related
%                       to the quadrature rule, or all of them (except boundary)
%                       as in the mesh structure from previous versions.
%     msh_evaluate_element_list: compute the parameterization (and its derivatives) at
%                       the quadrature points in a given list of elements.
%     msh_boundary_side_from_interior: create a msh object, with quadrature
%                       points only on one boundary side. Useful to apply boundary
%                       conditions in weak form.
%
% Copyright (C) 2009, 2010, 2014 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
% Copyright (C) 2014 Elena Bulgarello, Sara Frizziero
% Copyright (C) 2015 Rafael Vazquez
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

function msh = msh_cartesian (breaks, qn, qw, geo, varargin)

% For the 1D case
  if (~iscell (breaks))
    breaks = {breaks};
  end

  breaks = cellfun (@unique, breaks, 'UniformOutput', false);

  boundary = true;
  if (~isfield (geo, 'map_der2'))
    msh.der2 = false;
  else
    msh.der2 = true;
  end
  
  if (~isempty (varargin))
    if (~rem (length (varargin), 2) == 0)
      error ('msh_cartesian: options must be passed in the [option, value] format');
    end
    for ii=1:2:length(varargin)-1
      if (strcmpi (varargin {ii}, 'der2'))
        msh.der2 = varargin {ii+1};
      elseif (strcmpi (varargin {ii}, 'boundary'))
        boundary = varargin {ii+1};
      else
        error ('msh_cartesian: unknown option %s', varargin {ii});
      end
    end
  end

% To be consistent, in the 1D case we rewrite everything as cell-arrays
  if (~iscell (qn))
    qn = {qn};
    qw = {qw};
  end
  if (~iscell (breaks))
    breaks = {breaks};
  end
  
  
  msh.ndim = numel (qn);
  if (isfield (geo,'rdim'))
    msh.rdim = geo.rdim;
  else
    for idim = 1:msh.ndim
      brk1{idim} = breaks{idim}(1);
    end
    msh.rdim = size (feval (geo.map, brk1), 1);
  end

  msh.qn = qn;
  msh.qw = qw;

  % this only for precaution, breaks should already not have any repetitions
  msh.breaks = cellfun (@unique, breaks, 'UniformOutput', false);

  msh.nel_dir = cellfun (@numel, msh.breaks) - 1;
  msh.nel = prod (msh.nel_dir);
  
  msh.nqn_dir = cellfun (@(x) size(x,1), qn);
  msh.nqn  = prod (msh.nqn_dir);
  
  if (boundary && msh.ndim > 1)
    for iside = 1:2*msh.ndim
%%    ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind = [2 2 1 1] in 2D;
      ind = setdiff (1:msh.ndim, ceil(iside/2)); 

      if (isempty (qw))
        msh.boundary(iside) = msh_cartesian (msh.breaks(ind), msh.qn(ind), [], geo.boundary(iside), 'boundary', false);
      else
        msh.boundary(iside) = msh_cartesian (msh.breaks(ind), msh.qn(ind), msh.qw(ind), geo.boundary(iside), 'boundary', false);
      end
      msh.boundary(iside).rdim = msh.rdim;
    end
  elseif (boundary && msh.ndim == 1)
    msh.boundary(1).ndim = 0;
    msh.boundary(2).ndim = 0;
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
      warning ('msh_cartesian: a function to compute second order derivatives has not been provided')
    end
  end
  
  msh = class (msh, 'msh_cartesian');

  if (~isempty (msh.qw))
    for idim = 1:msh.ndim
      dim_length = msh.breaks{idim}(end) - msh.breaks{idim}(1);
      if (abs (sum (msh.qw{idim}(:)) - dim_length) > 1e-10)
        warning ('geopdes:check_quadrature', 'msh_cartesian: inconsistent quadrature formula')
      end
    end
  end
  
end
