% MSH_2D: constructor of the class for 2d tensor product meshes.
%
%     msh = msh_2d (breaks, qn, qw, geometry, opts);
%
% INPUTS:
%     
%     breaks:   breaks along each direction in parametric space (repetitions are ignored)
%     qn:       quadrature nodes along each direction in parametric space
%     qw:       quadrature weights along each direction in parametric space
%     geometry: structure representing the geometrical mapping
%     opts:     if opts == 'no boundary', the boundary terms are not computed
%   
% OUTPUT:
%
%     msh: class containing the following fields and methods
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
%     boundary      (1 x 4 struct-array)    it contains a one-dimensional 'msh' structure for each edge of the boundary 
%
%     METHOD NAME
%     msh_evaluate_col: computes the parameterization (and its derivatives) of
%                       the quadrature points in one column of the mesh, i.e.,
%                       fixing the element in the first parametric direction.
%     msh_eval_boundary_side: computes the parameterization in one boundary side
%                       of the domain.
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

function msh = msh_2d (breaks, qn, qw, geo, opts)

  if (nargin < 5) 
    opts = '';
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
  
  if (~strcmpi (opts, 'no boundary'))
    for iside = 1:4
      ind = mod (floor ((iside+1)/2), 2) + 1;  %ind = [2 2 1 1];
      msh.boundary(iside).breaks = msh.breaks{ind};
      msh.boundary(iside).nel = size (qn{ind},2);
      msh.boundary(iside).nqn = size (qn{ind},1);
    end
  else
    msh.boundary = [];    
  end

  msh.map = geo.map;
  msh.map_der = geo.map_der;
  msh = class (msh, 'msh_2d');

end
