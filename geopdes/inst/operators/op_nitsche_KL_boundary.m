% OP_NITSCHE_KL_BOUNDARY: compute matrix to impose zero rotation
%
%   N_mat = op_nitsche_KL_boundary (spu_param, spv_param, msh_side, msh_side_fi, E_coeff, nu_coeff, thickness, C_penalty)
%
% INPUT:
%     
%    spu_param:   space structure in the paramatric domain (see sp_vector/sp_evaluate_col_param)
%    spv_param:   space structure in the paramatric domain (see sp_vector/sp_evaluate_col_param)
%    msh_side:    mesh structure on the boundary
%    msh_side_fi: mesh side from interior, with information on the normal derivative
%    E_coeff:     Young modulus evaluated at quadrature points
%    nu_coeff:    Poisson ratio evaluated at quadrature points
%    thickness:   thickness of the shell
%    C_penalty:   penalization term
%   
% OUTPUT:
%
%     N_mat:       the computed matrix, to be added in the left hand-side
%
% Copyright (C) 2023, 2024 Giuliano Guarino, Rafael Vazquez
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

function mat = op_nitsche_KL_boundary (spu_param, spv_param, msh_side, msh_side_from_interior, E_bnd, nu_bnd, thickness, penalty_coeff)

nel = msh_side.nel; 
nqn = msh_side.nqn;
ndim = msh_side_from_interior.ndim;
rdim = msh_side.rdim;

qw_jacdet    = msh_side.quad_weights .* msh_side.jacdet;
normals      = msh_side.normal;
tangents     = msh_side.geo_map_jac;
% tangents_der = msh_side.geo_map_der2;

nsh_u = spu_param.nsh_max;
nsh_v = spv_param.nsh_max;
% Preallocating rows, columns and values for the interface
M_rows    = zeros(nsh_v*nsh_u*nel,1);
M_columns = zeros(nsh_v*nsh_u*nel,1);
M_values  = zeros(nsh_v*nsh_u*nel,1);
% R_rows    = zeros(nsh_v      *nel,1);
% R_values  = zeros(nsh_v      *nel,1);

M_ncounter = 0;
% R_ncounter = 0;

% Penalty coefficient
mu_r = penalty_coeff*E_bnd(1,1)*thickness^3;

for ise = 1:nel
    K = zeros(nsh_v,nsh_u);
    % R = zeros(nsh_v,1);
    for iqn = 1:nqn
        w = qw_jacdet(iqn,ise);
        % Tangent and derivative
        T_t  = tangents    (:,1,iqn,ise);
        % T_tt = tangents_der(:,1,1,iqn,ise);
        % Geometry
        Aa    = msh_side_from_interior.geo_map_jac (:,:,iqn,ise);
        Aa_b  = msh_side_from_interior.geo_map_der2(:,:,:,iqn,ise);
        % % Exact rotation
        % if ~isempty(rot_coeffs)
        %     Ta = Ta_L(:,iqn,ise);
        % end
        % Test and Trial shape functions derivatives
        u_a   = permute(reshape(spu_param.shape_function_gradients (:,:,  iqn,:,ise), [rdim,ndim,nsh_u]),[1,3,2    ]);
        u_ab  = permute(reshape(spu_param.shape_function_hessians  (:,:,:,iqn,:,ise), [rdim,ndim,ndim,nsh_u]),[1,4,2,3  ]);
        v_a   = permute(reshape(spv_param.shape_function_gradients (:,:,  iqn,:,ise), [rdim,ndim,nsh_v]),[1,3,2    ]);
        v_ab  = permute(reshape(spv_param.shape_function_hessians  (:,:,:,iqn,:,ise), [rdim,ndim,ndim,nsh_v]),[1,4,2,3  ]);
        % Material and normal
        Young     = E_bnd (iqn,ise); 
        Poisson   = nu_bnd(iqn,ise);
        [A_ort,B_ort,C_ort,D_ort] = isotropic_stiffness (Young, Poisson, thickness);
        normal    = normals(:,iqn,ise);

        % Eventually flipping the tangent
        if (normal.'*cross(T_t,cross(Aa(:,1),Aa(:,2)))<0)
            T_t  = -T_t;
            % T_tt = -T_tt;
        end
        % Rotation and forces and moment fluxes
%        [Tu, Fu, Mu] = op_KL_fluxes(u_a, u_ab, u_abc, Aa, Aa_b, Aa_bc, T_t, T_tt, A_ort, B_ort, C_ort, D_ort, C_ort_e, D_ort_e, normal);
        [Tu, Mu] = op_KL_rotation_fluxes(u_a, u_ab, Aa, Aa_b, T_t, D_ort, normal);
        [Tv, Mv] = op_KL_rotation_fluxes(v_a, v_ab, Aa, Aa_b, T_t, D_ort, normal);

        % Contribution of one Gauss point to the stiffness matrix
        K = K + w*mu_r*Tv'*Tu;
%        R = R + w*mu_r*Tv'*Ta;
        K = K - w*Mv'*Tu - w*Tv'*Mu;
%        R = R - w*Mv'*Ta;
    end
    % Assembling in the global matrix
    % rows and column dof list of the element in the global numeration
    rows_v    = spv_param.connectivity(:,ise);
    columns_v = spu_param.connectivity(:,ise);
    [rows_s,columns_s] = ndgrid(rows_v,columns_v);
    % Global rows and columns sparse vector
    M_rows   ((M_ncounter+1) : (M_ncounter+nsh_v*nsh_u) ) = rows_s(:);
    M_columns((M_ncounter+1) : (M_ncounter+nsh_v*nsh_u) ) = columns_s(:);
    M_values ((M_ncounter+1) : (M_ncounter+nsh_v*nsh_u) ) = K(:);
%    R_rows   ((R_ncounter+1) : (R_ncounter+nsh_v))        = rows_v(:);
%    R_values ((R_ncounter+1) : (R_ncounter+nsh_v))        = R(:);
    % Updating the ncounter
    M_ncounter = M_ncounter + nsh_v*nsh_u;
%    R_ncounter = R_ncounter + nsh_v;
end

mat = sparse(M_rows,M_columns,M_values,spv_param.ndof,spu_param.ndof);
% rhs = full(sparse(R_rows,ones(size(R_rows)),R_values,global_start_index(2),1));

end