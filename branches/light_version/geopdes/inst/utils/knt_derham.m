% KNT_DERHAM: construct the knot vectors for discrete B-spline spaces in the De Rham diagram.
%
%  USAGE:
%
%     [knots, degree] = knt_derham (knots_h1, degree_h1, output_space)
%
% INPUTS:
%
%     knots_h1:     knot vectors for the H^1-conforming space, cell array of size (1 x ndim)
%     degree_h1:    degree for the H^1 space, vector of size (1 x ndim)
%     output_space: string, with one of the following: 'H1', 'Hcurl', 'Hdiv', 'L2'
%
% OUTPUT
%
%    The output are the knot vectors and the degrees for the chosen output space
%    Their size depend on the space:
%
% For H^1 and L^2 (scalar spaces)
%
%     NAME               TYPE          SIZE       DESCRIPTION
%     knots              cell-array  (1 x ndim)
%     degree             vector      (1 x ndim)
%
% For H(curl) and H(div) (vectorial spaces) each component of the
% cell-array contains the same information for one component of the space.
% That is, the knot vector is a cell-array of cell-arrays.
%
%     NAME               TYPE          SIZE    
%     knots              cell-array  (1 x ndim)
%     degree             cell-array  (1 x ndim)
%
% Copyright (C) 2010, 2015 Rafael Vazquez
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


function [knots, degree] = knt_derham (knots_h1, degree_h1, output_space)

  if (nargin < 3)
    output_space = 'Hcurl';
  end
  
  if (nargout > 2)
    error (['Wrong number of output parameters. The use of knt_derham has changed since version 2.1.' ...
    'Please read the help, or the examples in GeoPDEs, to understand its usage'])
  end
  
  if (~iscell (knots_h1))
    knots_h1 = {knots_h1};
  end

  ndim = numel (knots_h1);
  if (ndim ~= numel (degree_h1))
    error ('knt_derham: wrong input parameters. The dimension of the knot vector and the degree do not match')
  end
  
  for ii = 1:ndim
    knots_l2{ii} = knots_h1{ii}(2:end-1);
  end

  if (strcmpi (output_space, 'H1'))
    knots = knots_h1;
    degree = degree_h1;
  elseif (strcmpi (output_space, 'L2'))
    knots = knots_l2;
    degree = degree_h1 - 1;
  elseif (strcmpi (output_space, 'Hcurl'))
    if (ndim < 2 || ndim > 3)
      error ('knt_derham: Wrong dimension to compute the H(curl) space')
    end
    knots = cell (ndim, 1);
    degree = cell (ndim, 1);
    for idim = 1:ndim
      knots{idim} = knots_h1;
      degree{idim} = degree_h1;
      knots{idim}{idim} = knots_l2{idim};
      degree{idim}(idim) = degree_h1(idim) - 1;
    end
  elseif (strcmpi (output_space, 'Hdiv'))
    if (ndim < 2 || ndim > 3)
      error ('knt_derham: Wrong dimension to compute the H(div) space')
    end
    knots = cell (ndim, 1);
    degree = cell (ndim, 1);
    for idim = 1:ndim
      knots{idim} = knots_l2;
      degree{idim} = degree_h1 - 1;
      knots{idim}{idim} = knots_h1{idim};
      degree{idim}(idim) = degree_h1(idim);
    end
  end
  
  
%   if (numel (knots_h1) == 2)
%     knots_u1 = {knots_l2{1} knots_h1{2}};
%     knots_u2 = {knots_h1{1} knots_l2{2}};
% 
%     if (nargout == 2)
%       varargout = {knots_u1 knots_u2};
%     elseif (nargout == 4 && nargin == 2)
%       degree1 = [degree(1) - 1, degree(2)];
%       degree2 = [degree(1), degree(2) - 1];
%       varargout = {knots_u1 knots_u2 degree1 degree2};
%     elseif (nargout == 3 && nargin == 1)
%       varargout = {knots_u1 knots_u2 knots_l2};
%     elseif (nargout == 6 && nargin == 2)
%       degree1 = [degree(1) - 1, degree(2)];
%       degree2 = [degree(1), degree(2) - 1];
%       degree_l2 = [degree(1) - 1, degree(2) - 1];
%       varargout = {knots_u1 knots_u2 knots_l2 degree1 degree2 degree_l2};
%     else
%       error ('knt_derham: wrong number of input or output parameters');
%     end
% 
%   elseif (numel (knots_h1) == 3)
%     knots_u1 = {knots_l2{1} knots_h1{2} knots_h1{3}};
%     knots_u2 = {knots_h1{1} knots_l2{2} knots_h1{3}};
%     knots_u3 = {knots_h1{1} knots_h1{2} knots_l2{3}};
% 
%     if (nargout == 3)
%       varargout = {knots_u1 knots_u2 knots_u3};
%     elseif (nargout == 6 && nargin == 2)
%       degree1 = [degree(1) - 1, degree(2), degree(3)];
%       degree2 = [degree(1), degree(2) - 1, degree(3)];
%       degree3 = [degree(1), degree(2), degree(3) - 1];
%       varargout = {knots_u1 knots_u2 knots_u3 degree1 degree2 degree3};
%     else
%       error ('knt_derham: wrong number of input or output parameters');
%     end
%   end
% 
end
