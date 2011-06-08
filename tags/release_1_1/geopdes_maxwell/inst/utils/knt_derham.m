% KNT_DERHAM: construct the knot vectors for discrete B-spline spaces in the De Rham diagram.
%
%  USAGE: 2D case
%
%     [knots1, knots2] = knt_derham (knots)
%     [knots1, knots2, degree1, degree2] = knt_derham (knots, degree)
%     [knots1, knots2, knots_mul] = knt_derham (knots)
%     [knots1, knots2, knots_mul, degree1, degree2, degree_mul] = knt_derham (knots, degree)
%
%  USAGE: 3D case
%
%     [knots1, knots2, knots3] = knt_derham (knots)
%     [knots1, knots2, knots3, degree1, degree2, degree3] = knt_derham (knots, degree)
%
% INPUTS:
%
%     knots:  knot vectors for the H^1-conforming space
%     degree: corresponding degree
%
% OUTPUT: 2D case
%
%     knots1:     knot vectors for the first component of the H(curl)-conforming space
%     degree1:    degree for the previous space
%     knots2:     knot vectors for the second component of the H(curl)-conforming space
%     degree2:    degree for the previous space
%     knots_mul:  knot vectors for the L^2-conforming discrete space
%     degree_mul: degree for the previous space
%
% OUTPUT: 3D case
%
%     knots1:     knot vectors for the first component of the H(curl)-conforming space
%     degree1:    degree for the previous space
%     knots2:     knot vectors for the second component of the H(curl)-conforming space
%     degree2:    degree for the previous space
%     knots3:     knot vectors for the third component of the H(curl)-conforming space
%     degree3:    degree for the previous space
%
% Copyright (C) 2010 Rafael Vazquez
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


function varargout = knt_derham (varargin)

  knots_h1 = varargin{1};
  if (nargin == 2)
    degree = varargin{2};
  end

  for ii = 1:numel(knots_h1)
    knots_l2{ii} = knots_h1{ii} (2:end-1);
  end

  if (numel (knots_h1) == 2)
    knots_u1 = {knots_l2{1} knots_h1{2}};
    knots_u2 = {knots_h1{1} knots_l2{2}};

    if (nargout == 2)
      varargout = {knots_u1 knots_u2};
    elseif (nargout == 4 && nargin == 2)
      degree1 = [degree(1) - 1, degree(2)];
      degree2 = [degree(1), degree(2) - 1];
      varargout = {knots_u1 knots_u2 degree1 degree2};
    elseif (nargout == 3 && nargin == 1)
      varargout = {knots_u1 knots_u2 knots_l2};
    elseif (nargout == 6 && nargin == 2)
      degree1 = [degree(1) - 1, degree(2)];
      degree2 = [degree(1), degree(2) - 1];
      degree_l2 = [degree(1) - 1, degree(2) - 1];
      varargout = {knots_u1 knots_u2 knots_l2 degree1 degree2 degree_l2};
    else
      error ('knt_derham: wrong number of input or output parameters');
    end

  elseif (numel (knots_h1) == 3)
    knots_u1 = {knots_l2{1} knots_h1{2} knots_h1{3}};
    knots_u2 = {knots_h1{1} knots_l2{2} knots_h1{3}};
    knots_u3 = {knots_h1{1} knots_h1{2} knots_l2{3}};

    if (nargout == 3)
      varargout = {knots_u1 knots_u2 knots_u3};
    elseif (nargout == 6 && nargin == 2)
      degree1 = [degree(1) - 1, degree(2), degree(3)];
      degree2 = [degree(1), degree(2) - 1, degree(3)];
      degree3 = [degree(1), degree(2), degree(3) - 1];
      varargout = {knots_u1 knots_u2 knots_u3 degree1 degree2 degree3};
    else
      error ('knt_derham: wrong number of input or output parameters');
    end
  end

end
