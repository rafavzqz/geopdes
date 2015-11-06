% SP_EVAL_STRESS_MSH: Evaluate the stress of a function, given by its degrees of freedom, at the points given by a msh object.
%
%   [eu, F] = sp_eval_stress_msh (u, space, msh, lambda, mu);
%
% INPUT:
%     
%     u:          vector of dof weights
%     space:      object defining the discrete space (see sp_bspline_3d)
%     msh:        object defining the points where to evaluate (see msh_3d)
%     lambda, mu: Lame' parameters
%
% OUTPUT:
%
%     eu: the gradient of the function evaluated in the given points 
%     F:  grid points in the physical domain, that is, the mapped points
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

function [stress, F] = sp_eval_stress_msh (u, space, msh, lambda, mu);

  if (space.ncomp == 1)
    error ('sp_eval_stress_msh: the stress is not computed for scalars')
  end

  ndim = numel (msh.qn);

  F  = zeros (ndim, msh.nqn, msh.nel);
  stress = zeros (space.ncomp, ndim, msh.nqn, msh.nel);

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    sp_col  = sp_evaluate_col (space, msh_col, 'value', false, ...
			       'gradient', true, 'divergence', true);

    uc_iel = zeros (size (sp_col.connectivity));
    uc_iel(sp_col.connectivity~=0) = ...
          u(sp_col.connectivity(sp_col.connectivity~=0));
    w_grad = repmat (reshape (uc_iel, [1, 1, 1, sp_col.nsh_max, msh_col.nel]), ...
                                  [sp_col.ncomp, ndim, msh_col.nqn, 1, 1]);
    w_div = repmat (reshape (uc_iel, [1, sp_col.nsh_max, msh_col.nel]), ...
                                  [msh_col.nqn, 1, 1]);

    sp_col.shape_function_gradients = ...
           reshape (sp_col.shape_function_gradients, sp_col.ncomp, ndim, ...
            msh_col.nqn, sp_col.nsh_max, msh_col.nel);
    sp_col.shape_function_divs = reshape (sp_col.shape_function_divs,  ...
                                      msh_col.nqn, sp_col.nsh_max, msh_col.nel);

    F(:,:,msh_col.elem_list) = msh_col.geo_map;

    for idim = 1:ndim
      x{idim} = reshape (msh_col.geo_map(idim,:,:), msh_col.nqn, msh_col.nel);
    end

    aux_grad = ...
       reshape (sum (w_grad .* sp_col.shape_function_gradients, 4), ...
                             sp_col.ncomp, ndim, msh_col.nqn, msh_col.nel);
    aux_div = reshape (sum (w_div .* sp_col.shape_function_divs, 2), ...
                             msh_col.nqn, msh_col.nel);

    mu_col = reshape (mu (x{:}), 1, 1, 1, msh_col.nel);
    stress(:,:,:,msh_col.elem_list) = bsxfun (@times, mu_col, ...
                         (aux_grad + permute (aux_grad, [2 1 3 4])));
    for idim = 1:ndim
      stress(idim,idim,:,msh_col.elem_list) = ...
        stress(idim,idim,:,msh_col.elem_list) + ...
        reshape (lambda (x{:}) .* aux_div, 1, 1, 1, msh_col.nel);
    end
  end

end
