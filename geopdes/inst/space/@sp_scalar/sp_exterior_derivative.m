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

  assert (strcmpi(space1.transform, 'grad-preserving'), ...
    'The first space cannot be the one for integral-preserving splines (or n-forms)')

  ndim = numel (space1.knots);
  if (ndim == 1)
    grad_curl = 'grad';
    assert (space1.degree == space2.degree+1, 'The degrees are not compatible')
    assert (numel(space1.knots{1}) == numel(space2.knots{1})+2, 'The knot vectors are not compatible')
  elseif (ndim == 2)
    if (strcmpi (space2.transform, 'curl-preserving'))
      grad_curl = 'grad';
      deg_shift = {[1 0], [0 1]};
      knt_shift = {[2 0], [0 2]};
    elseif (strcmpi (space2.transform, 'div-preserving'))
      grad_curl = 'curl';
      deg_shift = {[0 1], [1 0]};
      knt_shift = {[0 2], [2 0]};
    else
      error ('The second space should be either curl-preserving or div-preserving')
    end
    for idim = 1:ndim
      assert (all(space1.degree == space2.scalar_spaces{idim}.degree+deg_shift{idim}), 'The degrees are not compatible')
      assert (all(cellfun(@numel,space1.knots) == (cellfun(@numel, space2.scalar_spaces{idim}.knots)+knt_shift{idim})), ...
        'The knot vectors are not compatible')
    end
  elseif (ndim == 3)
    grad_curl = 'grad';
    deg_shift = {[1 0 0], [0 1 0], [0 0 1]};
    knt_shift = {[2 0 0], [0 2 0], [0 0 2]};
    assert (strcmpi(space2.transform, 'curl-preserving'), ...
      'The second space should be the one for curl-conforming splines')
    for idim = 1:ndim
      assert (all(space1.degree == space2.scalar_spaces{idim}.degree+deg_shift{idim}), 'The degrees are not compatible')
      assert (all(cellfun(@numel,space1.knots) == (cellfun(@numel, space2.scalar_spaces{idim}.knots)+knt_shift{idim})), ...
        'The knot vectors are not compatible')
    end
  end
  
  diff_ops = op_geom_exterior (space1.knots, space1.degree, grad_curl);
  diff_op = diff_ops{1};
end
