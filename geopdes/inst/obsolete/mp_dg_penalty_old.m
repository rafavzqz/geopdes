% MP_DG_PENALTY: apply the interior penalty method between patches. To be
%  used in multipatch geometries with RT and NDL spaces.
%
% USAGE:
%
%  A = mp_dg_penalty (space_v, msh, gnum, ornt, interfaces, visc, Cpen)
%
% INPUT:
%
%    space_v:     cell-array with space objects for the velocity (see sp_vector_div_transform)
%    msh:         cell-array with mesh object (see msh_cartesian)
%    gnum:        numbering of the degrees of freedom for each patch (see mp_interface_hdiv)
%    ornt:        orientation of the degrees of freedom (see mp_interface_hdiv)
%    interfaces:  interface information (see mp_geo_load)
%    visc:        function handle to compute the viscosity. For now it is
%                     assumed to be the same for all patches
%    Cpen:        penalization constant
%
% OUTPUT:
%
%    A: computed matrix, to be added to the global matrix (see mp_solve_stokes_div_conforming)
%
%   For more details, see:
%    J.A. Evans, T.J.R. Hughes
%    Isogeometric divergence-conforming B-splines for the Darcy-Stokes-Brinkman equations
%    Math. Models Meth. Appl. Sci., 2012
%
% Copyright (C) 2014, 2015 Rafael Vazquez
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

function A = mp_dg_penalty_old (space, msh, gnum, ornt, interfaces, visc, Cpen)

if (isa (space, 'sp_multipatch'))
  warning ('For spaces of the class SP_MULTIPATCH, using MP_DG_PENALTY')
  A = mp_dg_penalty (space, msh, interfaces, visc, Cpen);
  return
end

ndof = max([gnum{:}]);
rA = []; cA = []; vA = [];

ndim = msh{1}.ndim;
rdim = msh{1}.rdim;

for iref = 1:numel(interfaces)
  patch(1) = interfaces(iref).patch1;
  patch(2) = interfaces(iref).patch2;
  side(1) = interfaces(iref).side1;
  side(2) = interfaces(iref).side2;

  msh_side(1) = msh_eval_boundary_side (msh{patch(1)}, side(1));
  msh_side(2) = msh_eval_boundary_side (msh{patch(2)}, side(2));
  msh_side_int    = msh_boundary_side_from_interior (msh{patch(1)}, side(1));
  msh_side_int(2) = msh_boundary_side_from_interior (msh{patch(2)}, side(2));

  sp_aux = space{patch(1)}.constructor (msh_side_int(1));
  sp_bnd(1) = struct (sp_precompute (sp_aux, msh_side_int(1), 'value', true, 'gradient', true));
  sp_aux = space{patch(2)}.constructor (msh_side_int(2));
  sp_bnd(2) = struct (sp_precompute (sp_aux, msh_side_int(2), 'value', true, 'gradient', true));
  
  [sp_bnd(2), msh_side(2)] = reorder_elements_and_quad_points (sp_bnd(2), msh_side(2), interfaces(iref), ndim);
  
% I assume there are no dicontinuities in the viscosity
  for idim = 1:rdim
    x{idim} = reshape (msh_side(1).geo_map(idim,:,:), msh_side(1).nqn, msh_side(1).nel);
  end
  coeff_at_qnodes = visc (x{:});

  charlen = (sum (msh_side(1).charlen, 1).^(1/(ndim-1)));
  charlen = repmat (charlen, msh_side(1).nqn, 1);
  
  for ii = 1:2
    for jj = 1:2
      [rB, cB, vB] = op_gradu_v_otimes_n (sp_bnd(jj), sp_bnd(ii), msh_side(ii), coeff_at_qnodes/2);
      [rC, cC, vC] = op_u_otimes_n_v_otimes_n (sp_bnd(jj), sp_bnd(ii), msh_side(jj), msh_side(ii), ...
          coeff_at_qnodes ./ charlen);
      vB = ornt{patch(ii)}(rB)' .* vB .* ornt{patch(jj)}(cB)';
      vC = ornt{patch(ii)}(rC)' .* vC .* ornt{patch(jj)}(cC)' * Cpen;
      rB = gnum{patch(ii)}(rB); cB = gnum{patch(jj)}(cB);
      rC = gnum{patch(ii)}(rC); cC = gnum{patch(jj)}(cC);

      rA = [rA rB cB rC];
      cA = [cA cB rB cC];
      vA = [vA -vB -vB vC];
    end
  end
end

A = sparse (rA, cA, vA, ndof, ndof);
end


function [sp_bnd, msh_side] = reorder_elements_and_quad_points (sp_bnd, msh_side, interface, ndim)

  if (ndim == 2)
    if (interface.ornt == -1)
      elem = msh_side.nel:-1:1;
      qpoints = msh_side.nqn:-1:1;
    else
      elem = 1:msh_side.nel;
      qpoints = 1:msh_side.nqn;  
    end
  elseif (ndim == 3)
    elem = reshape (1:msh_side.nel, msh_side.nel_dir);
    qpoints = reshape (1:msh_side.nqn, msh_side.nqn_dir);
    if (interface.flag == -1)
      elem = elem';
      qpoints = qpoints';
    end
    if (interface.ornt1 == -1)
      elem = flipud (elem);
      qpoints = flipud (qpoints);
    end
    if (interface.ornt2 == -1)
      elem = fliplr (elem);
      qpoints = fliplr (qpoints);
    end
    elem = elem(:)';
    qpoints = qpoints(:)';
  end

  msh_side.charlen      = msh_side.charlen(qpoints,elem);
  msh_side.normal       = msh_side.normal(:,qpoints,elem);
  msh_side.jacdet       = msh_side.jacdet(qpoints,elem);
  msh_side.geo_map      = msh_side.geo_map(:,qpoints,elem);
  msh_side.geo_map_jac  = msh_side.geo_map_jac(:,:,qpoints,elem);
  msh_side.quad_weights = msh_side.quad_weights(qpoints,elem);
  msh_side.quad_nodes   = msh_side.quad_nodes(:,qpoints,elem);
      
  sp_bnd.nsh                      = sp_bnd.nsh(elem);
  sp_bnd.connectivity             = sp_bnd.connectivity(:,elem);
  sp_bnd.shape_functions          = sp_bnd.shape_functions(:,qpoints,:,elem);
%   sp_bnd.shape_function_divs      = sp_bnd.shape_function_divs(qpoints,:,elem);
  sp_bnd.shape_function_gradients = sp_bnd.shape_function_gradients(:,:,qpoints,:,elem);

end
