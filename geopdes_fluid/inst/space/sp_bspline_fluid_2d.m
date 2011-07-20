% SP_BSPLINE_FLUID_2D: Construct different pair of B-Splines spaces on the physical domain for fluid problems.
%
%   [spv, spp, PI] = sp_bspline_fluid_2d (elem_name, knots, nsub, ...
%                                               degreep, regularity, msh)
%
% INPUTS:
%
%   elem_name:  the name of the element. Right now 'TH' (Taylor-Hood), 
%                'NDL' (Nedelec, 2nd family), 'RT' (Raviart-Thomas) and
%                'SG' (SubGrid) are supported.
%   knots:      knot vector of the coarse geometry.
%   nsub:       number of subdivisions of each interval.
%   degreep:    degree of the pressure space along each parametric direction
%   regularity: continuity of the pressure space along each parametric direction
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

function [spv, spp, PI] = sp_bspline_fluid_2d (element_name, ...
                   knots, nsub_p, degree_p, regularity_p, msh)

% The pressure space is the same in the four cases
knotsp = kntrefine (knots, nsub_p-1, degree_p, regularity_p);
spp = sp_bspline_2d_phys (knotsp, degree_p, msh);

% Construction of the knot vectors and the discrete space for the velocity
switch (lower (element_name))
  case {'th'}
    degree_v = degree_p + 1;
    regularity_v = regularity_p;
    nsub_v = nsub_p;
    knots_v = kntrefine (knots, nsub_v-1, degree_v, regularity_v);

    spv = do_spv_component_wise_2d__ (knots_v, degree_v, knots_v, degree_v, msh);
    spv.spfun = @(MSH) do_spv_component_wise_2d__ (knots_v, degree_v, knots_v, degree_v, MSH);

    PI = speye (spp.ndof);
  case {'sg'}
    degree_v = degree_p + 1;
    regularity_v = regularity_p+1;
    nsub_v = 2*nsub_p;
    knots_v = kntrefine (knots, nsub_v-1, degree_v, regularity_v);

    spv = do_spv_component_wise_2d__ (knots_v, degree_v, knots_v, degree_v, msh);
    spv.spfun = @(MSH) do_spv_component_wise_2d__ (knots_v, degree_v, knots_v, degree_v, MSH);

    PI = speye (spp.ndof);
  case {'ndl'}
    degree_v1 = degree_p+1;
    degree_v2 = degree_p+1;
    knots_v1{1}  = [knotsp{1}(1) knotsp{1} knotsp{1}(end)];
    knots_v1{2}  = sort ([knotsp{2}, unique(knotsp{2})]);
    knots_v2{1}  = sort ([knotsp{1}, unique(knotsp{1})]);
    knots_v2{2}  = [knotsp{2}(1) knotsp{2} knotsp{2}(end)];

    spv = do_spv_piola_transform_2d__ (knots_v1, degree_v1, knots_v2, degree_v2, msh);
    spv.spfun = @(MSH) do_spv_piola_transform_2d__ (knots_v1, degree_v1, knots_v2, degree_v2, MSH);

    PI = speye (spp.ndof);
  case {'rt'}
    degree_v1 = [degree_p(1)+1 degree_p(2)];
    degree_v2 = [degree_p(1) degree_p(2)+1];
    knots_v1{1}  = [knotsp{1}(1) knotsp{1} knotsp{1}(end)];
    knots_v1{2}  = knotsp{2};
    knots_v2{1}  = knotsp{1};
    knots_v2{2}  = [knotsp{2}(1) knotsp{2} knotsp{2}(end)];

    spv = do_spv_piola_transform_2d__ (knots_v1, degree_v1, knots_v2, degree_v2, msh);
    spv.spfun = @(MSH) do_spv_piola_transform_2d__ (knots_v1, degree_v1, knots_v2, degree_v2, MSH);

    PI = b2nst__ (spp, knotsp, degree_p, msh);
  otherwise
    error ('sp_bspline_fluid_2d: unknown element type')
end

end

