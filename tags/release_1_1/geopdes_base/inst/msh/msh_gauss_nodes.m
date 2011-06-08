% MSH_GAUSS_NODES: define nodes and weights for a multi-dimensional
%                                 tensor-product Gauss quadrature rule.
%
%   rule = msh_gauss_nodes (nnodes)
%
% INPUT:
%
%     nnodes:   number of qadrature nodes along each direction
%
% OUTPUT:
%
%  rule:    cell array containing the nodes and weights of the rule  
%           along each direction (rule{idir}(1,:) are the nodes and 
%           rule{idir}(2,:) the weights)
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
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


function rule = msh_gauss_nodes (nnodes)
  
  for idir = 1:numel (nnodes)
    [rule{idir}(1,:), rule{idir}(2,:)] = grule (nnodes(idir));
  end

end
