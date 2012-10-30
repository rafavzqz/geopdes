% TSPLINE_BOUNDARY: construct the GeoPDEs mesh and space boundary structures from Bezier extraction.
%
%     [msh_bnd, space_bnd] = tspline_boundary (tspline, ibnd);
%
% INPUTS:
%
%     tspline: Bezier extraction read from a file of the T-splines plug-in for Rhino
%                (see read_bezier_extraction)
%     ibnd:    boundary side, as ordered in the "boundary" field of tspline
%   
% OUTPUT:
%
%     msh_bnd: structure containing the following fields
%
%        FIELD_NAME    (SIZE)                    DESCRIPTION
%        nel           (scalar)                  number of elements of the partition
%        nqn           (scalar)                  maximum number of quadrature nodes per element
%        nqn_elem      (1 x nel vector)          actual number of quadrature nodes on each element
%        quad_weights  (nqn x nel vector)        weights associated to the quadrature nodes
%        geo_map       (rdim x nqn x nel vector) physical coordinates of the quadrature nodes
%        geo_map_jac   
%               (rdim x ndim x nqn x nel vector) Jacobian matrix of the map evaluated at the quadrature nodes
%        jacdet        (nqn x nel)               determinant of the Jacobian evaluated at the quadrature points
%
%     space: structure containing the following fields
%
%        FIELD_NAME      (SIZE)                  DESCRIPTION
%        ndof            (scalar)                total number of degrees of freedom of the T-spline space (not only the boundary)
%        nsh_max         (scalar)                maximum number of shape functions per element
%        nsh             (1 x nel vector)        actual number of shape functions on each element
%        connectivity    (nsh_max x nel vector)  indices of basis functions that do not vanish on each element
%        shape_functions (nqn x nsh_max x nel)   basis functions evaluated at each quadrature node on each element
%
%  Warning: the space structure uses the global numbering, and the field dofs is not computed.
%    This is different from the other examples in GeoPDEs.
%
% Copyright (C) 2012 Rafael Vazquez
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

function [msh_bnd, space_bnd] = tspline_boundary (tspline, ibnd)

% Construct the boundary fields (only for "side" sets)
  if (tspline.ndim == 3 || tspline.rdim == 3)
    error ('Until now the boundary is only implemented for planar surfaces')
  elseif (~strcmp (tspline.boundary{ibnd}.type, 'side'))
    error ('The boundary is only computed for side sets')
  end

% Check whether the degree is the same for all the Bezier elements or not
  boundary_elements = tspline.elements(tspline.boundary{ibnd}.elem);
  check_degree = @(x) all (x == boundary_elements(1).degree);
  same_degree = all (cellfun (check_degree, {boundary_elements.degree}));

  if (~same_degree)
    error ('The case of different degrees on each Bezier element has not been implemented yet')
  end

  nsides = tspline.boundary{ibnd}.nsides;

% Initialize to zero, to allocate memory.
%  The values assigned to nqn and nsh_max may be changed (updated) later
  msh_bnd.nel = nsides;
  msh_bnd.nqn = prod (boundary_elements(1).degree(1:end-1) + 1);
  msh_bnd.nqn_elem = zeros (1, nsides);
  msh_bnd.quad_weights = zeros (msh_bnd.nqn, nsides);
  msh_bnd.geo_map = zeros (tspline.rdim, msh_bnd.nqn, nsides);
  msh_bnd.geo_map_jac = zeros (tspline.rdim, tspline.ndim, msh_bnd.nqn, nsides);
  msh_bnd.jacdet = zeros (msh_bnd.nqn, nsides);

  space_bnd.ncomp = 1;
  space_bnd.ndof = tspline.ndof;
  space_bnd.nsh_max = prod (boundary_elements(1).degree(1:end-1) + 1);
  space_bnd.nsh = zeros (1, msh_bnd.nel);
  space_bnd.connectivity = zeros (space_bnd.nsh_max, msh_bnd.nel);
  space_bnd.shape_functions = zeros (msh_bnd.nqn, space_bnd.nsh_max, msh_bnd.nel);

% Precompute the quadrature points and univariate Bernstein polynomials 
%  in the parent element  
  degree = boundary_elements(1).degree;
  nqn_dir = boundary_elements(1).degree + 1;

  rule = msh_gauss_nodes (nqn_dir);
  for idir = 1:tspline.ndim
    [qn{idir}, qw{idir}] = deal (rule{idir}(1,:)', rule{idir}(2,:)');
  end

  for idir = 1:tspline.ndim
    knt = [-ones(1, degree(idir)+1), ones(1, degree(idir)+1)];
    s = findspan (degree(idir)+1, degree(idir), qn{idir}, knt);
    bernstein_univ = basisfunder (s, degree(idir), qn{idir}, knt, 1);
    shape_funs{idir} = reshape (bernstein_univ(:,1,:), nqn_dir(idir), degree(idir)+1);
    shape_fun_ders{idir} = reshape (bernstein_univ(:,2,:), nqn_dir(idir), degree(idir)+1);
  end    

  for iside = 1:nsides
    element = boundary_elements(iside);
    degree = element.degree;
% Identify the univariate Bernstein polynomial from the position,
%   L(eft), R(ight), B(ottom) and T(op).
    switch tspline.boundary{ibnd}.position{iside}
      case {'L'}
        bernstein_indices = sub2ind (degree+1, ones (1, degree(2)+1), 1:degree(2)+1);
        ind = 2;
      case {'R'}
        bernstein_indices = sub2ind (degree+1, (degree(1)+1)*ones (1, degree(2)+1), 1:degree(2)+1);
        ind = 2;
      case {'B'}
        bernstein_indices = sub2ind (degree+1, 1:degree(1)+1, ones (1, degree(1)+1));
        ind = 1;
      case {'T'}
        bernstein_indices = sub2ind (degree+1, 1:degree(1)+1, (degree(2)+1)*ones (1, degree(1)+1));
        ind = 1;
    end

% Compute the reduced extraction operator
    C = element.extraction;
    bsp_indices = any (C(:,bernstein_indices), 2);
    C_bnd = C(bsp_indices, bernstein_indices);
    connectivity_iside = element.connectivity(bsp_indices);
    nsh_iside = numel (connectivity_iside);

    ctrl_points = tspline.control_points(1:tspline.rdim, connectivity_iside);
    weights = tspline.control_points(4, connectivity_iside);

    W = sparse (diag (weights));

    bernstein_shape_funs = shape_funs{ind};
    bernstein_shape_fun_der = shape_fun_ders{ind};
% Univariate non-rational basis functions and derivatives on the boundary
    bsp_shape_funs = C_bnd * bernstein_shape_funs';
    bsp_shape_fun_tang_der =  C_bnd * bernstein_shape_fun_der';
% Univariate rational basis functions and derivatives on the boundary
    Wb = weights * bsp_shape_funs;
    nrb_shape_funs = bsxfun (@rdivide, W * bsp_shape_funs, Wb);
    dWb = weights * bsp_shape_fun_tang_der;
    nrb_shape_fun_tang_der = ...
      bsxfun (@rdivide, W * bsp_shape_fun_tang_der, Wb) - ...
      bsxfun (@times, W * bsp_shape_funs, dWb ./ (Wb.^2));

% Apply the parameterization: quadrature points in the physical domain
    geo_map = ctrl_points * nrb_shape_funs;
    geo_map_jac = ctrl_points * nrb_shape_fun_tang_der;
    jacdet = geopdes_norm__ (geo_map_jac);

% Write the fields into the mesh and space structures
    nqn_iside = numel (qw{ind});

    msh_bnd.nqn_elem(iside) = nqn_iside;
    msh_bnd.quad_weights(1:nqn_iside,iside) = qw{ind};
    msh_bnd.geo_map(:,1:nqn_iside,iside) = geo_map;
    msh_bnd.jacdet(1:nqn_iside,iside) = jacdet;

    space_bnd.nsh(iside) = nsh_iside;
    space_bnd.connectivity(1:nsh_iside,iside) = connectivity_iside';
    space_bnd.shape_functions(1:nqn_iside,1:nsh_iside,iside) = permute (nrb_shape_funs, [2 1]);
   
  end
  
% Set the number of quadrature points, and the number of basis functions
%  per element, as the maximum
  msh_bnd.nqn = max (msh_bnd.nqn_elem);
  space_bnd.nsh_max = max (space_bnd.nsh);
  
end
