% msh_partition: partititon a mesh in parts suitable for parallel execution
%
%  [breaks boundaries patch] = mp_msh_partition( pieces, data ...)
%
% INPUT:
%
% pieces:   number of required parts
% data:	    a list of data in the [data, value] format. Possible data are
%               'breaks', partition trying to divide the elements evenly
%               'knots', if also 'degree' is present partition trying to distribute 
%                        degrees of freedom evenly, else fallback to distributing elements
%               'degree', gives the degree of the spline functions
%               'index', if set return onl the breaks for the part having that label else return a
%                        cell array with all the parts.
%
% OUTPUT:
%
% breaks:     a cell array containing the breaks along all dimension of all the parts
% boundaries: a cell array containing for each part the logical vector of
%             needed boundaries
% patch:      which patch do we work in
%
% Copyright (C) 2011 Andrea Bressan
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 2 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.


function [breaks boundaries patch] = mp_msh_partition(pieces, knots, degree, piece_id)

ptch_num = numel(knots);
if (pieces< ptch_num)
  error('there must be at least a processor per patch');
end

patch      = 1 + mod(piece_id -1, ptch_num);
sub_id     = ceil(piece_id/ptch_num);
sub_pieces = floor (pieces/ptch_num);
if patch <= mod(pieces, ptch_num)
  sub_pieces=sub_pieces+1;
end
[breaks boundaries] = msh_partition(sub_pieces, 'knots', knots{patch}, 'degree', degree, 'index', sub_id);
end