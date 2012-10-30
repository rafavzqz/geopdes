% TSPLINE_MESH_SPACE: construct the GeoPDEs mesh and space structures from Bezier extraction.
%
%     [msh, space] = tspline_mesh_space (tspline);
%
% INPUTS:
%     
%     tspline: Bezier extraction read from a file of the T-splines plug-in for Rhino
%                (see read_bezier_extraction)
%   
% OUTPUT:
%
%     msh: structure containing the following fields
%
%        FIELD_NAME    (SIZE)                    DESCRIPTION
%        nel           (scalar)                  number of elements of the partition
%        nqn           (scalar)                  number of quadrature nodes per element
%        quad_weights  (nqn x nel vector)        weights associated to the quadrature nodes
%        geo_map       (rdim x nqn x nel vector) physical coordinates of the quadrature nodes
%        geo_map_jac   (rdim x ndim x nqn x nel) Jacobian matrix of the map evaluated at the quadrature nodes
%        jacdet        (nqn x nel)               determinant of the Jacobian evaluated at the quadrature points
%
%     space: structure containing the following fields
%
%        FIELD_NAME      (SIZE)                        DESCRIPTION
%        ndof            (scalar)                      total number of degrees of freedom
%        nsh_max         (scalar)                      maximum number of shape functions per element
%        nsh             (1 x nel vector)              actual number of shape functions on each element  
%        connectivity    (nsh_max x nel vector)        indices of basis functions that do not vanish on each element
%        shape_functions (nqn x nsh_max x nel)         basis functions evaluated at each quadrature node on each element
%        shape_function_gradients
%                        (ndim x nqn x nsh_max x nel)  basis function gradients evaluated at each quadrature node on each element
%
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2012 Rafael Vazquez
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

function [msh, space] = tspline_mesh_space (tspline, option, npts)

  if (nargin < 2)
    option = 'quadrature';
  end
  if (nargin < 3 && strcmpi (option, 'quadrature'))
    npts = [];
  elseif (nargin < 3 && strcmpi (option, 'plot'))
    npts = repmat (10, 1, tspline.ndim);
  end

  msh.nel = tspline.nel;
  space.ndof = tspline.ndof;
  space.ncomp = 1;

  connectivity = cellfun (@transpose, {tspline.elements.connectivity}, ...
      'UniformOutput', false);
  space.nsh = cellfun (@numel, connectivity);
  space.nsh_max = max (space.nsh);

% Initialize nqn and initialize all the arrays to zero, to allocate memory
  all_degrees = unique ([tspline.elements.degree]);
  if (isempty (npts))
    npts = (max (all_degrees) + 1) * ones (1, tspline.ndim);
  end
  msh.nqn = prod (npts);

  if (strcmpi (option, 'quadrature'))
    msh.quad_weights = zeros (msh.nqn, msh.nel);
  end
  msh.geo_map = zeros (tspline.rdim, msh.nqn, msh.nel);
  msh.geo_map_jac = zeros (tspline.rdim, tspline.ndim, msh.nqn, msh.nel);
  msh.jacdet = zeros (msh.nqn, msh.nel);
  
  space.connectivity = zeros (space.nsh_max, msh.nel);
  space.shape_functions = zeros (msh.nqn, space.nsh_max, msh.nel);
  space.shape_function_gradients = zeros (tspline.ndim, msh.nqn, space.nsh_max, msh.nel);

% Precompute the Bernstein polynomials at the quadrature points
%  for all the degrees that are present in the T-spline space.
%  We compute everything in the parent element [-1 1]^n.
  if (strcmpi (option, 'quadrature'))
    rule = msh_gauss_nodes (npts);
    for idir = 1:tspline.ndim
      [qn{idir}, qw{idir}] = deal (rule{idir}(1,:)', rule{idir}(2,:)');
    end
    if (tspline.ndim == 2)
      quad_weights = kron (qw{2}, qw{1});
      msh.quad_weights = repmat (quad_weights(:), 1, msh.nel);
    elseif (tspline.ndim == 3)
      quad_weights = kron (qw{3}, kron (qw{2}, qw{1}));
      msh.quad_weights = repmat (quad_weights(:), 1, msh.nel);
    end
  elseif (strcmpi (option, 'plot'))
    for idir = 1:tspline.ndim
      qn{idir} = linspace (-1+eps, 1-eps, npts(idir));
    end        
  end

% Univariate Bernstein polynomials, computed as B-spline functions
  for ideg = all_degrees
    for idir = 1:tspline.ndim
      knt = [-ones(1, ideg+1), ones(1, ideg+1)];
      s = ideg * ones (size (qn{idir})); % s = findspan (ideg, ideg, qn{idir}, knt);
      bernstein_univ = basisfunder (s, ideg, qn{idir}, knt, 1);
      shape_funs{ideg}{idir} = reshape (bernstein_univ(:,1,:), npts(idir), ideg + 1);
      shape_fun_ders{ideg}{idir} = reshape (bernstein_univ(:,2,:), npts(idir), ideg + 1);
    end
  end

% This is a weird way to decide which degrees combinations are present in the space
aux = {tspline.elements.degree};
mult_deg = cellfun (@str2num, unique (cellfun (@num2str, aux, 'UniformOutput', false)), ...
        'UniformOutput', false);
  
% Multivariate Bernstein polynomials
  if (tspline.ndim == 2)
    for ii = 1:numel(mult_deg)
      ideg = mult_deg{ii}(1); jdeg = mult_deg{ii}(2);
      bernstein_shape_funs{ideg, jdeg} = kron (shape_funs{jdeg}{2}, shape_funs{ideg}{1});
      bernstein_shape_fun_grad{ideg, jdeg}{1} = kron (shape_funs{jdeg}{2}, shape_fun_ders{ideg}{1});
      bernstein_shape_fun_grad{ideg, jdeg}{2} = kron (shape_fun_ders{jdeg}{2}, shape_funs{ideg}{1});
    end
  elseif (tspline.ndim == 3)
    for ii = 1:numel(mult_deg)
      ideg = mult_deg{ii}(1); jdeg = mult_deg{ii}(2); kdeg = mult_deg{ii}(3);
      bernstein_shape_funs = kron (shape_funs{kdeg}{3}, kron (shape_funs{jdeg}{2}, shape_funs{ideg}{1}));
      bernstein_shape_fun_grad{ideg, jdeg}{1} = kron (shape_funs{kdeg}{3}, kron (shape_funs{jdeg}{2}, shape_fun_ders{ideg}{1}));
      bernstein_shape_fun_grad{ideg, jdeg}{2} = kron (shape_funs{kdeg}{3}, kron (shape_fun_ders{jdeg}{2}, shape_funs{ideg}{1}));
      bernstein_shape_fun_grad{ideg, jdeg}{3} = kron (shape_fun_ders{kdeg}{3}, kron (shape_funs{jdeg}{2}, shape_funs{ideg}{1}));
    end
  end

% For every Bezier element, we start computing the basis functions in the
%   parent element [-1 1]^n
    for iel = 1:msh.nel
      element = tspline.elements(iel);
      deg = element.degree;
      C = element.extraction;
      ctrl_points = tspline.control_points(1:tspline.rdim, element.connectivity);
      weights = tspline.control_points(4, element.connectivity);

      W = sparse (diag (weights));

      bernstein_shape_funs = kron (shape_funs{deg(2)}{2}, shape_funs{deg(1)}{1});
      bernstein_shape_fun_grad{1} = kron (shape_funs{deg(2)}{2}, shape_fun_ders{deg(1)}{1});
      bernstein_shape_fun_grad{2} = kron (shape_fun_ders{deg(2)}{2}, shape_funs{deg(1)}{1});
      
% Non-rational T-spline basis functions and derivatives in the parent element
      bsp_shape_funs = C * bernstein_shape_funs';
      for idim = 1:tspline.ndim
        bsp_shape_fun_grad{idim} = C *bernstein_shape_fun_grad{idim}';
      end
% Rational T-spline basis functions and derivatives in the parent element
      Wb = weights * bsp_shape_funs;
      nrb_shape_funs = bsxfun (@rdivide, W * bsp_shape_funs, Wb);
      for idim = 1:tspline.ndim
        dWb = weights * bsp_shape_fun_grad{idim};
        nrb_shape_fun_grad{idim} = ...
          bsxfun (@rdivide, W * bsp_shape_fun_grad{idim}, Wb) - ...
          bsxfun (@times, W * bsp_shape_funs, dWb ./ (Wb.^2));
      end

% Apply the parameterization: quadrature points in the physical domain
      geo_map = ctrl_points * nrb_shape_funs;
      geo_map_jac = zeros (tspline.ndim, tspline.rdim, msh.nqn);
      for idim = 1:tspline.ndim
        geo_map_jac(:,idim,:) = ctrl_points * nrb_shape_fun_grad{idim};
      end
      if (tspline.ndim == tspline.rdim)
        jacdet = geopdes_det__ (geo_map_jac);
      else
        error ('3D surfaces (or other manifolds) are not implemented yet')
      end

% Write the fields into the msh structure
      msh.geo_map(:,:,iel) = geo_map;
      msh.geo_map_jac(:,:,:,iel) = geo_map_jac;
      msh.jacdet(:,iel) = abs (jacdet);

% Write the fields into the space structure
%  The basis functions and derivatives are still in the parent element
      space.connectivity(1:space.nsh(iel),iel) = element.connectivity';
      space.shape_functions(:,1:space.nsh(iel),iel) = permute (nrb_shape_funs, [2 1]);
      for idim = 1:tspline.ndim
        space.shape_function_gradients(idim,:,1:space.nsh(iel),iel) = ...
            permute (nrb_shape_fun_grad{idim}, [2 1]);
      end
    end %iel

% Compute the derivatives in the physical domain, applying the transformation
    if (tspline.ndim == tspline.rdim)
      JinvT = geopdes_invT__ (msh.geo_map_jac);
      JinvT = reshape (JinvT, [tspline.ndim, tspline.ndim, msh.nqn, msh.nel]);
      space.shape_function_gradients = ...
          geopdes_prod__ (JinvT, space.shape_function_gradients);
    else
      error ('3D surfaces (or other manifolds) are not implemented yet')
    end

end
