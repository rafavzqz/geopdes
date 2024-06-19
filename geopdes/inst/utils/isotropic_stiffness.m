function [A_ort,B_ort,C_ort,D_ort] = isotropic_stiffness(young,poisson,thickness)
%
% Copyright (C) 2024 Giuliano Guarino
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


% Computing the stiffness matrix in two dimensions ------------------------
E = young/(1-poisson^2)*[1,poisson,0;poisson,1,0;0,0,0.5-poisson/2];    
% membrane and flexural stiffness matrices --------------------------------
A_ort = E*thickness;
D_ort = E*thickness*thickness*thickness/12;
% transforming to tensor --------------------------------------------------
A_ort = stiffness_mat2ten(A_ort);
D_ort = stiffness_mat2ten(D_ort);
B_ort = zeros(2,2,2,2);
C_ort = zeros(2,2,2,2);
end




function T = stiffness_mat2ten(M)
% this function transform a plate stiffness matrix in a tensor 
    T(1,1 , 1,1) = M(1,1);
    T(1,1 , 1,2) = M(1,3);
    T(1,1 , 2,1) = M(1,3);
    T(1,1 , 2,2) = M(1,2);

    T(1,2 , 1,1) = M(3,1);
    T(1,2 , 1,2) = M(3,3);
    T(1,2 , 2,1) = M(3,3);
    T(1,2 , 2,2) = M(3,2);

    T(2,1 , 1,1) = M(3,1);
    T(2,1 , 1,2) = M(3,3);
    T(2,1 , 2,1) = M(3,3);
    T(2,1 , 2,2) = M(3,2);

    T(2,2 , 1,1) = M(2,1);
    T(2,2 , 1,2) = M(2,3);
    T(2,2 , 2,1) = M(2,3);
    T(2,2 , 2,2) = M(2,2);
end