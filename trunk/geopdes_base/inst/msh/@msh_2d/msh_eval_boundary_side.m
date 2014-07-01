% MSH_EVAL_BOUNDARY_SIDE: evaluate the parameterization in one boundary side of the domain.
%
%     msh_side = msh_eval_boundary_side (msh, iside);
%
% INPUTS:
%     
%    msh:   mesh object (see msh_2d)
%    iside: number of the boundary side to compute, from 1 to 4 (see the file geo_specs for documentation about edge numbering)
%   
% OUTPUT:
%
%     msh_side: structure that contains the following fields
%
%     FIELD_NAME    (SIZE)                  DESCRIPTION
%     side_number   (scalar)                number of the side
%     nel           (scalar)                number of elements of the boundary side
%     nqn           (scalar)                number of quadrature nodes per element
%     quad_nodes    (2 x nqn x nel vector)  coordinates of the quadrature nodes in parametric domain
%     quad_weights  (nqn x nel vector)      weights associated to the quadrature nodes
%     geo_map       (2 x nqn x nel vector)  physical coordinates of the quadrature nodes
%     geo_map_jac   (2 x 2 x nqn x nel)     Jacobian matrix of the map evaluated at the quadrature nodes
%     jacdet        (nqn x nel)             determinant of the Jacobian evaluated at the quadrature points
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2014 Rafael Vazquez
% Copyright (C) 2014 Adriano Cortes
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

function msh_side = msh_eval_boundary_side (msh, iside)

  ind = mod (floor ((iside+1)/2), 2) + 1;  %ind = [2 2 1 1];
  ind2 = floor ((iside+1)/2);              %ind2 = [1 1 2 2];

  msh_side = msh.boundary(iside);

  msh_side.quad_weights = msh.qw{ind};
  msh_side.quad_nodes = zeros ([2, size(msh.qn{ind})]);
  msh_side.quad_nodes(ind,:,:) = msh.qn{ind};
  if (iside == 2 || iside == 4)
    msh_side.quad_nodes(ind2,:,:) = 1;
  end

  qn1 = msh_side.quad_nodes(1,:,:);
  qn2 = msh_side.quad_nodes(2,:,:);
  F   = feval (msh.map, [qn1(:), qn2(:)]');
  jac = feval (msh.map_der, [qn1(:), qn2(:)]');

  msh_side.geo_map = reshape (F, size (msh_side.quad_nodes));
  msh_side.geo_map_jac = reshape(jac, 2, 2, msh_side.nqn, msh_side.nel);
  msh_side.jacdet = geopdes_norm__ (squeeze (msh_side.geo_map_jac(:,ind,:,:)));
  clear F jac
  
  [JinvT, ~] = geopdes_invT__ (msh_side.geo_map_jac);
  JinvT = reshape (JinvT, [2, 2, msh_side.nqn, msh_side.nel]);
  normal = zeros (2, msh_side.nqn, msh_side.nel);
  normal(ind2,:,:) = (-1)^iside;
  normal = reshape (normal, [2, msh_side.nqn, 1, msh_side.nel]);
  normal = geopdes_prod__ (JinvT, normal);
  normal = reshape (normal, [2, msh_side.nqn, msh_side.nel]);
  norms = repmat (reshape (geopdes_norm__ (normal), [1, msh_side.nqn, msh_side.nel]), [2 1 1]);
  %Now normalize
  normal = normal ./ norms;
  msh_side.normal = normal;
  clear JinvT;

% Computation of the normal characteristic length for Nitsche's method
% Computed as in "Weak Dirichlet boundary condition for wall-bounded
%  turbulent flows", Y. Bazilevs and T.J.R. Hughes (2007).
  if(iside == 1)
      xi_span_charlen = 0.5*(msh.breaks{1}(2) - msh.breaks{1}(1)).*ones(1,msh.nel_dir(2));
      eta_span_charlen = 0.5*diff(msh.breaks{2});
  elseif(iside == 2)
      xi_span_charlen = 0.5*(msh.breaks{1}(end) - msh.breaks{1}(end-1)).*ones(1,msh.nel_dir(2));
      eta_span_charlen = 0.5*diff(msh.breaks{2});
  elseif(iside == 3)
      xi_span_charlen = 0.5*diff(msh.breaks{1});
      eta_span_charlen = 0.5*(msh.breaks{2}(2) - msh.breaks{2}(1)).*ones(1,msh.nel_dir(1));
  else
      xi_span_charlen = 0.5*diff(msh.breaks{1});
      eta_span_charlen = 0.5*(msh.breaks{2}(end) - msh.breaks{2}(end-1)).*ones(1,msh.nel_dir(1));
  end

  geo_map_jac = msh_side.geo_map_jac;
  
  for iel = 1:msh_side.nel
      for nquad = 1:msh_side.nqn
          geo_map_jac(:,1,nquad,iel) = xi_span_charlen(iel).*geo_map_jac(:,1,nquad,iel);
          geo_map_jac(:,2,nquad,iel) = eta_span_charlen(iel).*geo_map_jac(:,2,nquad,iel);
      end
  end
  
  [Jinv, ~] = geopdes_inv__ (geo_map_jac);
  Jinv = reshape (Jinv, [2, 2, msh_side.nqn, msh_side.nel]);
  normal = reshape (normal, [2, msh_side.nqn, 1, msh_side.nel]);
  normal = geopdes_prod__ (Jinv, normal);
  normal = reshape (normal, [2, msh_side.nqn, msh_side.nel]);
  norms = geopdes_norm__ (normal);
  msh_side.charlen = (2./norms);
  clear normal norms

end
