% MSH_SET_QUAD_NODES: compute location and weights of quadrature nodes in a tensor product grid.
%
%  [qn, qw] = msh_set_quad_nodes (breaks, rule)
%  [qn, qw] = msh_set_quad_nodes (breaks, rule, limits)
%
% INPUT:
%     
%  breaks:  breaks along each direction (in parametric space)
%  rule:    cell array containing the nodes and weights of the rule to be 
%           used in each direction (rule{idir}(1,:) are the nodes and rule{idir}(2,:) the weights)
%  limits:  (optional) limits of the interval in which the quadrature rule 
%            was defined. By default it is [-1 1].
%
% OUTPUT:
%
%   qn: {qnu, [qnv, [qnw]]}  quadrature nodes
%   qw: {qwu, [qwv, [qww]]}  quadrature weights
% 
% Copyright (C) 2009, 2010 Carlo de Falco
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

function [qn, qw] = msh_set_quad_nodes (knots, rule, limits)

  if (nargin < 3)
    limits = [-1 1];
  end

  ndir = numel (rule);
  for idir = 1:ndir
    uniqueknots  = unique(knots{idir});
    uniqueknotsl = uniqueknots(1:end-1);
    du           = diff (uniqueknots);
  
    nel = length (uniqueknots)-1;

    [p1, w1] = deal (rule{idir}(1,:)', rule{idir}(2,:)');
    nqn = numel (p1);

    qn{idir} = repmat (uniqueknotsl, nqn, 1) + repmat (p1-limits(1), 1, nel) .* repmat (du/diff(limits), nqn, 1);
    qw{idir} = repmat (w1, 1, nel) .* repmat (du/diff(limits), nqn, 1);
  end

end

