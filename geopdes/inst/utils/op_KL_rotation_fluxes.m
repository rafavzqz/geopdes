% OP_KL_ROTATION_FLUXES: compute the rotation and the bending moment along the boundary identified by normal and tangent
%
%   [Rotation, Moment] = op_KL_rotation_fluxes(u_a, u_ab, Aa, Aa_b, T_t, Dabcd_ort, normal)
%
% INPUT:
%
%  u_a:        [3 x nsh x 2] multidimenional array. u_a(i,j,k) is the derivative along xi_k of the j-th shape function associated to the i-th Cartesian component of the dispacement field
%  u_ab:       [3 x nsh x 2 x 2]  multidimenional array. u_a(i,j,k,l) is the derivative along xi_l of u_a(i,j,k)
%  Aa:         [3 x 2] matrix.  A_a(i,j) is the derivative along xi_j of the i-th Cartesian component of the surface map
%  Aa_b:       [3 x 2 x 2] multidimenional array. A_a(i,j,k) is the derivative along xi_k of A_a(i,j)
%  T_t:        [3 x 1] vector. T_t(i) is the derivative of the i-th Cartesian component of the map of the boundary with respect to the curvilinear cordinate that parametrizes the boundary
%  Dabcd_ort:  [2 x 2 x 2 x 2] multidimensional array. Constitutive tensor in local orthonormal coordinates (normalized first covariant and second contravariant basis vectors). 
%  normal:     [3 x 1] vector. Outer unit normal 
%  
%
% OUTPUT:
%
%  Rotation:   [1 x nsh] vector. Multiplied for the vector of the correspondent dofs returns the rotation along the tangent
%  Moment:     [1 x nsh] vector. Multiplied for the vector of the correspondent dofs returns the moment along the tangent
%
%
% NOTES:
%
%  1) While normal has to be normalized (unit vector) T_t does not have to
%  2) Dabcd_ort is in tensor form. If the constitutive relationship is available in Voigt form it has to be converted
%  3) The quantities involved all refer to a single point along the boundary line
%
%
% Copyright (C) 2024 Giuliano Guarino
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.
%
%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.



function [Rotation, Moment] = op_KL_rotation_fluxes(u_a, u_ab, Aa, Aa_b, T_t, Dabcd_ort, normal)
% Data from input =========================================================
% Number of shape functions -----------------------------------------------
nsh = size(u_a,2);
% =========================================================================



%% Geometry ===============================================================
% Covariant base on the shell surface -------------------------------------
A1 = Aa(:,1);
A2 = Aa(:,2);
% Vector V orthogonal to the surface (not normalized) and its derivative --
V   = cross(A1,A2);
% Norm of V and its derivatives -------------------------------------------
lam      = norm(V);
% Normalized vector orthogonal to the surface and its derivatives ---------
A3        = V/lam;
% Carthesian basis --------------------------------------------------------
e1_cart = [1;0;0];
e2_cart = [0;1;0];
e3_cart = [0;0;1];
% Auxiliary matrix for linearization of n0 and derivatives ----------------
auxM1 = [cross(e1_cart,A1),cross(e2_cart,A1),cross(e3_cart,A1)];
auxM2 = [cross(e1_cart,A2),cross(e2_cart,A2),cross(e3_cart,A2)];
auxM  = A3*A3'/lam;
% Covariant components of the surface metric tensor -----------------------
A11_cov = A1'*A1;
A12_cov = A1'*A2;
A21_cov = A2'*A1;
A22_cov = A2'*A2;
Aab_cov = [A11_cov,A12_cov;A21_cov,A22_cov];
% Contravariant components of the surface metric tensor -------------------
Aab_con = Aab_cov^-1;
% Contravariant basis on the surface ---------------------------------------
Aa_con(:,1) = Aab_con(1,1)*Aa(:,1) + Aab_con(1,2)*Aa(:,2);
Aa_con(:,2) = Aab_con(2,1)*Aa(:,1) + Aab_con(2,2)*Aa(:,2);
% Vector t and its derivatives --------------------------------------------
tau = norm(T_t);
t   = T_t/tau;
% Checking the orientation of the normal vector ---------------------------
aux = normal'*cross(t,A3);
if aux>0
    flip = false;
elseif aux<0
    flip = true;
end
% Vector n and its derivatives --------------------------------------------
if flip
    n    = -( cross(t,A3) );
else
    n    = cross(t,A3);
end
% Cov components of the vector n ------------------------------------------
n1 = n'*A1;
n2 = n'*A2;
na = [n1,n2];
% Cross products useful for the computations ------------------------------
Aa_bxAc   = zeros(3,2,2,2);
for i1 = 1:2
    for i2 = 1:2
        for i3 = 1:2
            Aa_bxAc(:,i1,i2,i3) = cross(Aa_b(:,i1,i2),Aa(:,i3));
        end
    end
end
% Orthonormal basis vectors and their derivatives -------------------------
e1 = A1/norm(A1);
e2 = cross(A3,e1);
ea = [e1,e2];
% =========================================================================



%% Displacements and rotations ============================================
% Rotation ----------------------------------------------------------------
txA3 = cross(t,A3);
txA3xA1 = cross(txA3,A1);
txA3xA2 = cross(txA3,A2);
phin = txA3xA1'/lam*u_a(:,:,2) - txA3xA2'/lam*u_a(:,:,1) + txA3'*auxM*(-auxM2*u_a(:,:,1)+auxM1*u_a(:,:,2));
% =========================================================================



%% Deformations ===========================================================
kab   = zeros(nsh,2,2);
% Cycle for computing eab and kab -----------------------------------------
for i1 = 1:2
    for i2 = 1:2
        kab(:,i1,i2) = -A3'*u_ab(:,:,i1,i2) - Aa_bxAc(:,i1,i2,1)'/lam*u_a(:,:,2) + Aa_bxAc(:,i1,i2,2)'/lam*u_a(:,:,1)+...
                       -Aa_b(:,i1,i2)'*auxM*auxM1*u_a(:,:,2) + Aa_b(:,i1,i2)'*auxM*auxM2*u_a(:,:,1);
    end
end
% =========================================================================



%% Stiffness coefficient ==================================================
% Computation of the cosine directors and derivatives ---------------------
eaAb_con = zeros(2,2);
for i1 = 1:2
    for i2 = 1:2
        eaAb_con(i1,i2) = ea(:,i1)'*Aa_con(:,i2);
    end
end
% Computation of the stiffness coefficients -------------------------------
Dabcd = zeros(2,2,2,2);
for i1 = 1:2 
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                for i5 = 1:2
                    for i6 = 1:2
                        for i7 = 1:2
                            for i8 = 1:2
                                Dabcd(i1,i2,i3,i4) = Dabcd(i1,i2,i3,i4) + Dabcd_ort(i5,i6,i7,i8)*eaAb_con(i5,i1)*eaAb_con(i6,i2)*eaAb_con(i7,i3)*eaAb_con(i8,i4);
                            end
                        end
                    end
                end
            end
        end
    end
end
% =========================================================================



%% Stress resultant =======================================================
mab    = zeros(nsh,2,2);
% Cycle for computing mab -------------------------------------------------
for i1 = 1:2
    for i2 = 1:2
        for i3 = 1:2
            for i4 = 1:2
                mab(:,i1,i2) = mab(:,i1,i2) + Dabcd(i1,i2,i3,i4)*kab(:,i3,i4);
            end
        end
    end
end
% =========================================================================



%% Forces for the penalization ============================================
% Cycle for computing Mt --------------------------------------------------
Mt = zeros(nsh,1);
for i1 = 1:2
   for i2 = 1:2
       Mt = Mt + mab(:,i1,i2)*na(i1)*na(i2);
   end
end
Mt = Mt';
% Eventual flipping of mt -------------------------------------------------
if flip
    Mt = -Mt;
end
% =========================================================================



%% Output =================================================================
Rotation = phin;
Moment   = Mt;
% =========================================================================
end











