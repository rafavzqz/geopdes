% OP_GRADFU_GRADV_MP: assemble the matrix A = [a(i,j)], a(i,j) = (grad (f(w) u_j), grad v_i) = 
%  = (f(w) grad u_j, grad v_i) + (u_j grad (f(w)), grad v_i), 
%  for C^1 multipatch splines, with "f" a given function
%  (derivative also needed), and "w" a discrete solution defined on the
%  same space. Useful to compute nonlinear terms.
%
%   [mat1, mat2] = op_gradfu_gradv_mp (space, msh, uhat, f, df, [patches]);
%
% INPUT:
%
%   space:   object representing the space of trial functions (see sp_multipatch_C1)
%   msh:     object defining the domain partition and the quadrature rule (see msh_multipatch)
%   uhat:    degrees of freedom of "w", of size space.ndof x 1.
%   f:       function handle to compute f(w).
%   df:      function handle to compute df(w), the gradient of f.
%   patches: list of patches where the integrals have to be computed. By default, all patches are selected.
%
% OUTPUT:
%
%   mat1:   assembled matrix for (f(w) grad u_j, grad v_i)
%   mat2:   assembled matrix for (u_j grad (f(w)), grad v_i)
% 
% Copyright (C) 2023 Michele Torre, Rafael Vazquez
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

function [A, B] = op_gradfu_gradv_mp (space, msh, uhat, f, df, patch_list)

  if (nargin < 6)
    patch_list = 1:msh.npatch;
  end

  if ((space.npatch ~= msh.npatch))
    error ('op_gradfu_gradv_mp: the number of patches does not coincide')
  end

  A = spalloc (space.ndof, space.ndof, 3*space.ndof);
  B = spalloc (space.ndof, space.ndof, 6*space.ndof);

  for iptc = patch_list
    [Cpatch_u, Cpatch_cols_u] = sp_compute_Cpatch (space, iptc);
    u_ptc = Cpatch_u * uhat(Cpatch_cols_u);

    [Ap, Bp] = op_gradfu_gradv_tp (space.sp_patch{iptc}, msh.msh_patch{iptc}, u_ptc, f, df);

    A(Cpatch_cols_u,Cpatch_cols_u) = ...
          A(Cpatch_cols_u,Cpatch_cols_u) + Cpatch_u.' * Ap * Cpatch_u;

    B(Cpatch_cols_u,Cpatch_cols_u) = ...
          B(Cpatch_cols_u,Cpatch_cols_u) + Cpatch_u.' * Bp * Cpatch_u;
  end

end
