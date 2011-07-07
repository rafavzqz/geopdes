% SP_BSPLINE_FLUID_3D_PHYS: Construct different pair of B-Splines spaces on the physical domain for fluid problems.
%
%   [spv, spp, PI] = sp_bspline_fluid_3d_phys (knotsv1, degreev1, ...
%          knotsv2, degreev2, knotsv3, degreev3, knotsp, degreep, ...
%          msh, fun_transform, press_proj)
%
% INPUTS:
%
%   elem_name: the name of the element. Right now 'TH' (Taylor-Hood), 
%              and 'SG' (SubGrid) are supported.
%   knotsv1:   knot vector of the space for the first parametric component of 
%               the velocity along each parametric direction
%   degreev1:  degree of the space for the first parametric component of
%               the velocity along each parametric direction
%   knotsv2:   knot vector of the space for the second parametric component of 
%               the velocity along each parametric direction
%   degreev2:  degree of the space for the second parametric component of
%               the velocity along each parametric direction
%   knotsv3:   knot vector of the space for the third parametric component of 
%               the velocity along each parametric direction
%   degreev3:  degree of the space for the third parametric component of
%               the velocity along each parametric direction
%   knotsp:    knot vector of the pressure space along each parametric direction
%   degreep:   degree of the pressure space along each parametric direction
%   msh:       msh class containing (in the field msh.qn) the points 
%                along each parametric direction in the parametric 
%                domain at which to evaluate, i.e. quadrature points 
%                or points for visualization (see msh_3d)
%
% OUTPUT:
%
%   spv: object representing the discrete velocity function space (see sp_vector_3d)
%   spp: object representing the discrete pressure function space (see sp_bspline_3d)
%   PI:  a projection matrix for the application of boundary conditions
%         for Raviart-Thomas spaces. The identity matrix in all other cases.
%
%   For more details, see:
%      A.Buffa, C.de Falco, G. Sangalli, 
%      IsoGeometric Analysis: Stable elements for the 2D Stokes equation
%      IJNMF, 2010
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011 Andrea Bressan, Rafael Vazquez
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

function [spv, spp, PI] = sp_bspline_fluid_3d_phys (element_name, ...
              knotsv1, degreev1, knotsv2, degreev2, knotsv3, degreev3, ...
              knotsp, degreep, msh)

spp = sp_bspline_3d (knotsp, degreep, msh);

switch (lower (element_name))
  case {'th', 'sg'}
    sp1 = sp_bspline_3d (knotsv1, degreev1, msh);
    sp2 = sp_bspline_3d (knotsv2, degreev2, msh);
    sp3 = sp_bspline_3d (knotsv3, degreev3, msh);
    spv = sp_vector_3d (sp1, sp2, sp3, msh);

    PI = speye (spp.ndof);
  case {'ndl', 'rt'}
    error ('NDL and RT elements have not been implemented in 3D yet')
  otherwise
    error ('sp_bspline_fluid_3d_phys: unknown element type')
end

end
