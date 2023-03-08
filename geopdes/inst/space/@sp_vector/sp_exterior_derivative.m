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

  assert (~strcmpi(space1.transform, 'integral-preserving'), ...
    'The first space cannot be the one for integral-preserving splines (or n-forms)')

  ndim = numel (space1.scalar_spaces{1}.knots);
  if (ndim == 1)
    error ('Not implemented. For dimension 1, use the scalar spaces')
  elseif (ndim == 2)
    if (strcmpi (space1.transform, 'curl-preserving'))
      grad_curl = 'grad';
      output_der = 2;
      deg_shift = {[0 1], [1 0]};
      knt_shift = {[0 2], [2 0]};
      degree = space1.scalar_spaces{1}.degree + [1 0];
      knots = {space1.scalar_spaces{2}.knots{1}, space1.scalar_spaces{1}.knots{2}};
    elseif (strcmpi (space1.transform, 'div-preserving'))
      grad_curl = 'curl';
      output_der = 2;
      deg_shift = {[1 0], [0 1]};
      knt_shift = {[2 0], [0 2]};
      degree = space1.scalar_spaces{1}.degree + [0 1];
      knots = {space1.scalar_spaces{1}.knots{1}, space1.scalar_spaces{2}.knots{2}};
    else
      error ('The second space should be either curl-preserving or div-preserving')
    end
    assert (isa(space2, 'sp_scalar'), 'The two spaces are not compatible')
    for idim = 1:ndim
      assert (all(space1.scalar_spaces{idim}.degree == space2.degree+deg_shift{idim}), 'The degrees are not compatible')
      assert (all(cellfun(@numel,space1.scalar_spaces{idim}.knots) == (cellfun(@numel, space2.knots)+knt_shift{idim})), ...
        'The knot vectors are not compatible')
    end
  elseif (ndim == 3)
    grad_curl = 'grad';
    if (strcmpi (space1.transform, 'curl-preserving'))
      output_der = 2;
      deg_shift = {[-1 1 1], [1 -1 1], [1 1 -1]};
      knt_shift = {[-2 2 2], [2 -2 2], [2 2 -2]};
      degree = space1.scalar_spaces{1}.degree + [1 0 0];
      knots = {space1.scalar_spaces{2}.knots{1}, space1.scalar_spaces{1}.knots{2}, space1.scalar_spaces{1}.knots{3}};
      assert (isa(space2, 'sp_vector'), 'The two spaces are not compatible')
      for idim = 1:ndim
        assert (all(space1.scalar_spaces{idim}.degree == space2.scalar_spaces{idim}.degree+deg_shift{idim}), 'The degrees are not compatible')
        assert (all(cellfun(@numel,space1.scalar_spaces{idim}.knots) == (cellfun(@numel, space2.scalar_spaces{idim}.knots)+knt_shift{idim})), ...
          'The knot vectors are not compatible')
      end
    elseif (strcmpi (space1.transform, 'div-preserving'))
      output_der = 3;
      deg_shift = {[1 0 0], [0 1 0], [0 0 1]};
      knt_shift = {[2 0 0], [0 2 0], [0 0 2]};
      degree = space1.scalar_spaces{1}.degree + [0 1 1];
      knots = {space1.scalar_spaces{1}.knots{1}, space1.scalar_spaces{2}.knots{2}, space1.scalar_spaces{3}.knots{3}};
      assert (isa(space2, 'sp_scalar'), 'The two spaces are not compatible')
      for idim = 1:ndim
        assert (all(space1.scalar_spaces{idim}.degree == space2.degree+deg_shift{idim}), 'The degrees are not compatible')
        assert (all(cellfun(@numel,space1.scalar_spaces{idim}.knots) == (cellfun(@numel, space2.knots)+knt_shift{idim})), ...
          'The knot vectors are not compatible')
      end
    else
      error ('Only implemented for curl-preserving and div-preserving transforms')
    end
  end
  
  diff_ops = op_geom_exterior (knots, degree, grad_curl);
  diff_op = diff_ops{output_der};
end
