% PIOLA_TRANSFORM_GRAD__: compute the gradients for the div-conforming Piola transform.
%
%     shape_fun_grads = piola_transform_grad__ (space, msh, shp, shg)
%
% INPUT:
%     
%    space: structure representing the discrete function space (see sp_vector_div_transform/sp_evaluate_col).
%    msh: structure containing the quadrature information (see msh_cartesian/msh_evaluate_col)
%    shp: basis functions evaluated in the parametric domain
%    shg: gradient of the basis functions, in the parametric domain
%
% OUTPUT:
%
%    shape_fun_grads: gradient of the basis functions in the physical domain
%
% Copyright (C) 2015 Rafael Vazquez
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

function shape_fun_grads = piola_transform_grad__ (sp, msh, shp, shg)

  shape_fun_grads = zeros (sp.ncomp, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);

  DF  = reshape (msh.geo_map_jac, msh.rdim, msh.ndim, 1, msh.nqn, 1, msh.nel);
  DF2 = reshape (msh.geo_map_der2, msh.rdim, msh.ndim, msh.ndim, msh.nqn, 1, msh.nel);
  shp = reshape (shp, 1, sp.ncomp_param, 1, msh.nqn, sp.nsh_max, msh.nel);
  shg = reshape (shg, 1, sp.ncomp_param, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);

  Jinv = geopdes_inv__ (msh.geo_map_jac);
  Jinv  = reshape (Jinv, [1, msh.ndim, msh.rdim, msh.nqn, 1, msh.nel]);
  JinvT = permute (Jinv, [3 2 1 4 5 6]);
  jacdet = reshape (geopdes_det__ (msh.geo_map_jac), msh.nqn, 1, msh.nel);


% Derivatives of |DF|, divided by |DF| to simplify
  Ddet = zeros (msh.rdim, msh.nqn, 1, msh.nel);
  DF2_JinvT = sum (sum (bsxfun (@times, DF2, JinvT), 1), 2);
  JinvT = reshape (JinvT, [1, msh.rdim, msh.ndim, msh.nqn, 1, msh.nel]);
  for idim = 1:msh.rdim
    aux = JinvT(:,idim,:,:,:,:);
    Ddet(idim,:,:,:) = sum (bsxfun (@times, aux, DF2_JinvT), 3);
  end
  clear aux DF2_JinvT

  for icomp = 1:sp.ncomp
    aux2 = sum (bsxfun (@times, DF(icomp,:,:,:,:,:), shg), 2);
    aux3 = reshape (sum (bsxfun (@times, DF(icomp,:,:,:,:,:), shp), 2), [1, msh.nqn, sp.nsh_max, msh.nel]);

    for jdim = 1:msh.rdim
      aux = sum (bsxfun (@times, DF2(icomp,:,:,:,:,:), JinvT(:,jdim,:,:,:,:)), 3);
      term1 = reshape (sum (bsxfun (@times, aux, shp), 2), [msh.nqn, sp.nsh_max, msh.nel]);

      term2 = reshape (sum (bsxfun (@times, aux2, JinvT(:,jdim,:,:,:,:)), 3), [msh.nqn, sp.nsh_max, msh.nel]);
      
      term3 = reshape (bsxfun (@times, -Ddet(jdim,:,:,:), aux3), [msh.nqn, sp.nsh_max, msh.nel]);

      shape_fun_grads(icomp,jdim,:,:,:) = bsxfun (@rdivide, term1 + term2 + term3, jacdet);
    end
  end

end
