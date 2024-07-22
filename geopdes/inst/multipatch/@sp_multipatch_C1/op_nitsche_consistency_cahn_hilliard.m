% OP_NITSCHE_CONSISTENCY_CAHN_HILLIARD: assemble the matrix A = [a(i,j)], a(i,j) = (epsilon (grad v n)_j, Delta u_i), 
%  with n the normal vector to the boundary, using the same space for trial and test functions.
%
%   mat = op_nitsche_consistency_cahn_hilliard (space, msh, sides, [coeff]);
%
% INPUT:
%
%  space:  object representing the space of trial functions (see sp_multipatch_C1)
%  msh:    object defining the mesh (see msh_multipatch)
%  sides:  boundary sides on which to compute the integrals
%  coeff:  function handle to compute the epsilon coefficient. If empty, it is taken equal to one.
%
% OUTPUT:
%
%  mat:    assembled matrix
% 
% Copyright (C) 2023, 2024 Michele Torre, Rafael Vazquez
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

function A = op_nitsche_consistency_cahn_hilliard (space, msh, nmnn_sides, coeff)

  if (nargin == 3 || isempty(coeff))
    coeff = @(varargin) ones (size(varargin{1}));
  end

  if (~isempty(nmnn_sides))

    A =  spalloc (space.ndof, space.ndof, 3*space.ndof);

    for iref = nmnn_sides
      for bnd_side = 1:msh.boundaries(iref).nsides  
        iptc = msh.boundaries(iref).patches(bnd_side);
        iside = msh.boundaries(iref).faces(bnd_side);

        msh_side = msh_eval_boundary_side (msh.msh_patch{iptc}, iside ) ;
        msh_side_int = msh_boundary_side_from_interior (msh.msh_patch{iptc}, iside ) ;

        sp_side = space.sp_patch{iptc}.constructor (msh_side_int);
        sp_side = sp_precompute ( sp_side , msh_side_int , 'gradient' , true, 'laplacian', true );

        for idim = 1:msh.rdim
          x{idim} = reshape (msh_side.geo_map(idim,:,:), msh_side.nqn, msh_side.nel);
        end
        coe_side = coeff (x{:});
        tmp = op_gradv_n_laplaceu(sp_side ,sp_side ,msh_side, coe_side);

        [Cpatch, Cpatch_cols] = sp_compute_Cpatch (space, iptc);
        A(Cpatch_cols,Cpatch_cols) = ...
              A(Cpatch_cols,Cpatch_cols) + Cpatch.' * tmp * Cpatch;
      end
    end

  else
    A = [];
  end
end
