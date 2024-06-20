% SP_EXTERIOR_DERIVAITVE: computes the exterior derivative as a matrix with size 
%  given by the dimension of two consecutive spaces in the De Rham sequence.
%
%   diff_op = sp_exterior_derivative (space1, space2);
%
% INPUT:
%
%   space1:  domain space of the exterior derivative (number of columns)
%   space2:  image space of the exterior derivative (number of rows)
%
% OUTPUT:
%
%   diff_op: sparse matrix representation of the differential operator.
% 
% Copyright (C) 2020-2023 Bernard Kapidani, Rafael Vazquez
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

function diff_op = sp_exterior_derivative (space1, space2)

  assert (space1.npatch == space2.npatch, 'The spaces correspond to different domains')
  npatch = space1.npatch;

  diff_op = sparse (space2.ndof, space1.ndof);
  
  for iptc = 1:npatch
    Dmat = sp_exterior_derivative (space1.sp_patch{iptc}, space2.sp_patch{iptc});
    
    ndof1 = space1.sp_patch{iptc}.ndof;
    ndof2 = space2.sp_patch{iptc}.ndof;
    if (isa (space1.sp_patch{iptc}, 'sp_vector'))
      dofs_ornt = space1.dofs_ornt{iptc};
      Dornt_1 = spdiags (dofs_ornt(:), 0, ndof1, ndof1);
    else
      Dornt_1 = speye (ndof1, ndof1);
    end
    if (isa (space2.sp_patch{iptc}, 'sp_vector'))
      dofs_ornt = space2.dofs_ornt{iptc};
      Dornt_2 = spdiags (dofs_ornt(:), 0, ndof2, ndof2);
    else
      Dornt_2 = speye (ndof2, ndof2);
    end
    gnum1 = space1.gnum{iptc}(:);
    gnum2 = space2.gnum{iptc}(:);
    diff_op(gnum2, gnum1) = Dornt_2 * Dmat * Dornt_1;
  end

end
