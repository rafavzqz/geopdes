% MSH_3D: constructor of the class for 3d tensor product meshes.
%
%     msh = msh_3d (breaks, qn, qw, geometry, opts);
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
%     nel_dir       (1 x 3 vector)          number of elements in each parametric direction
%     nelcol        (scalar)                number of elements in one "column" of the mesh (actually, nel_dir(2)*nel_dir(3))
%     nqn           (scalar)                number of quadrature nodes per element
%     nqnu          (scalar)                number of quadrature nodes per element in the first parametric direction
%     nqnv          (scalar)                number of quadrature nodes per element in the second parametric direction
%     nqnw          (scalar)                number of quadrature nodes per element in the third parametric direction
%     breaks        (1 x 3 cell-array)      unique(breaks)
%     qn            (1 x 3 cell-array)      quadrature nodes along each direction in parametric domain
%     qw            (1 x 3 cell-array)      quadrature weights along each direction in parametric space
%     boundary      (1 x 6 struct-array)    it contains a one-dimensional 'msh' structure for each face of the boundary 
%
%     METHOD NAME
%     msh_evaluate_col: computes the parameterization (and its derivatives) of
%                       the quadrature points in one column of the mesh, i.e.,
%                       fixing the element in the first parametric direction.
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

function msh = msh_3d (breaks, qn, qw, geo, opts)

  if (nargin < 5) 
    opts = '';
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
  
  if (~strcmpi (opts, 'no boundary'))
    for iside = 1:6
%%    ind =[2 3; 2 3; 1 3; 1 3; 1 2; 1 2]
      ind = setdiff (1:3, ceil(iside/2)); 

      boundary.qn = {qn{ind}};
      boundary.qw = {qw{ind}};

      boundary.nel_dir = msh.nel_dir(ind);
      boundary.nel = prod (boundary.nel_dir);
      boundary.nqn_dir = msh.nqn_dir(ind);
      boundary.nqn = prod (boundary.nqn_dir);

      msh.boundary(iside) = boundary;
    end
  else
    msh.boundary = [];    
  end

  msh.map = geo.map;
  msh.map_der = geo.map_der;
  msh = class (msh, 'msh_3d');

end
