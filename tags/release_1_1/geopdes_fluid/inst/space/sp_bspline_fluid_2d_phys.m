% SP_BSPLINE_FLUID_2D_PHYS: Construct different pair of B-Splines spaces on the physical domain for fluid problems.
%
%   [spv, spp, PI] = sp_bspline_fluid_2d_phys (knotsv1, degreev1, ...
%          knotsv2, degreev2, knotsp, degreep, msh)
%
% INPUTS:
%
%   elem_name: the name of the element. Right now 'TH' (Taylor-Hood), 
%               'NDL' (Nedelec, 2nd family), 'RT' (Raviart-Thomas) and
%               'SG' (SubGrid) are supported.
%   knotsv1:   knot vector of the space for the first parametric component of 
%               the velocity along each parametric direction
%   degreev1:  degree of the space for the first parametric component of
%               the velocity along each parametric direction
%   knotsv2:   knot vector of the space for the second parametric component of 
%               the velocity along each parametric direction
%   degreev2:  degree of the space for the second parametric component of
%               the velocity along each parametric direction
%   knotsp:    knot vector of the pressure space along each parametric direction
%   degreep:   degree of the pressure space along each parametric direction
%   msh:       structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%
% OUTPUT:
%
%   spv: structure representing the discrete velocity function space
%   spp: structure representing the discrete pressure function space
%         see sp_bspline_2d_phys for more details
%   PI:  a projection matrix for the application of boundary conditions
%         for Raviart-Thomas spaces
%
%   For more details, see:
%      A.Buffa, C.de Falco, G. Sangalli, 
%      IsoGeometric Analysis: Stable elements for the 2D Stokes equation
%      IJNMF, 2010
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011 Rafael Vazquez
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

function [spv, spp, PI] = sp_bspline_fluid_2d_phys (element_name, ...
              knotsv1, degreev1, knotsv2, degreev2, knotsp, degreep, msh)

spp = sp_bspline_2d_phys (knotsp, degreep, msh, 'gradient', false);

switch (lower (element_name))
  case {'th', 'sg'}
    spv = do_spv_component_wise_2d__ (knotsv1, degreev1, knotsv2, degreev2, msh);
    spv.spfun = @(MSH) do_spv_component_wise_2d__ (knotsv1, degreev1, knotsv2, degreev2, MSH);
    PI = speye (spp.ndof);

  case {'ndl'}
    spv = do_spv_piola_transform_2d__ (knotsv1, degreev1, knotsv2, degreev2, msh);
    spv.spfun = @(MSH) do_spv_piola_transform_2d__ (knotsv1, degreev1, knotsv2, degreev2, MSH);
    PI = speye (spp.ndof);

  case {'rt'}
    spv = do_spv_piola_transform_2d__ (knotsv1, degreev1, knotsv2, degreev2, msh);
    spv.spfun = @(MSH) do_spv_piola_transform_2d__ (knotsv1, degreev1, knotsv2, degreev2, MSH);
    PI = b2nst__ (spp, knotsp, degreep, msh);

  otherwise
    error ('sp_bspline_fluid_2d_phys: unknown element type')
end

end

