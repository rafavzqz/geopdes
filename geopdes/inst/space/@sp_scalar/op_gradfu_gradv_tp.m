% OP_GRADFU_GRADV_TP: assemble the matrix A = [a(i,j)], a(i,j) = (grad (f(w) u_j), grad v_i) = 
%  = (f(w) grad u_j, grad v_i) + (u_j grad (f(w)), grad v_i), 
%  exploiting the tensor product structure, with "f" a given function
%  (derivative also needed), and "w" a discrete solution defined on the
%  same space. Useful to compute nonlinear terms.
%
%   [mat1, mat2] = op_gradfu_gradv_tp (space, msh, uhat, f, df);
%
% INPUT:
%
%   space:   object representing the space of trial and test functions (see sp_scalar)
%   msh:     object defining the domain partition and the quadrature rule (see msh_cartesian)
%   uhat:    degrees of freedom of "w", of size space.ndof x 1.
%   f:       function handle to compute f(w).
%   df:      function handle to compute df(w), the gradient of f.
%
% OUTPUT:
%
%   mat1:   assembled matrix for (f(w) grad u_j, grad v_i)
%   mat2:   assembled matrix for (u_j grad (f(w)), grad v_i)
%
% Copyright (C) 2011, Carlo de Falco, Rafael Vazquez
% Copyright (C) 2016, 2017, Rafael Vazquez
% Copyright (C) 2023, Michele Torre, Rafael Vazquez
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

function [A, B] = op_gradfu_gradv_tp (space, msh, uhat, f, df)

  for idim = 1:msh.ndim
    size1 = size (space.sp_univ(idim).connectivity);
    if (size1(2) ~= msh.nel_dir(idim))
      error ('The discrete space is not associated to the mesh')
    end
  end

  A = spalloc (space.ndof, space.ndof, 3*space.ndof);
  B = spalloc (space.ndof, space.ndof, 6*space.ndof);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'gradient', true);

    % Evaluate the field and its gradient at the Gaussian points
    utemp = sp_eval_msh (uhat, sp_col, msh_col, {'value', 'gradient'});
    u = utemp{1};
    gradu = utemp{2};

    % Polynomial formulation for the double-well
    coeffs_A = f(u);
    A = A + op_gradu_gradv (sp_col, sp_col, msh_col, coeffs_A);

    coeffs_B = df(u);
    coeffs_Bv = gradu;
    for idim = 1:msh.ndim
      coeffs_Bv(idim,:,:) = coeffs_Bv(idim,:,:) .* reshape(coeffs_B, 1, size(coeffs_B,1), size(coeffs_B,2));
    end
    B = B + op_vel_dot_gradu_v (sp_col, sp_col, msh_col, coeffs_Bv).';
  end
