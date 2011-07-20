% SP_BSPLINE_FLUID_2D_PHYS: Construct different pair of B-Splines spaces on the physical domain for fluid problems.
%
%   [spv, spp, PI] = sp_bspline_fluid_2d_phys (elem_name, knots, nsub, ...
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
%   msh:        msh object containing (in the field msh.qn) the points 
%                along each parametric direction in the parametric 
%                domain at which to evaluate, i.e. quadrature points 
%                or points for visualization (see msh_2d).
%
% OUTPUT:
%
%   spv: object representing the discrete velocity function space (see sp_vector_2d, sp_vector_2d_piola_transform)
%   spp: object representing the discrete pressure function space (see sp_bspline_2d)
%   PI:  a projection matrix for the application of boundary conditions
%         for Raviart-Thomas spaces. The identity matrix in all other cases.
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
                   knots, nsub_p, degree_p, regularity_p, msh)

% The pressure space is the same in the four cases
knotsp = kntrefine (knots, nsub_p-1, degree_p, regularity_p);
spp = sp_bspline_2d (knotsp, degree_p, msh);

% Construction of the knot vectors for the velocity components
switch (lower (element_name))
  case {'th'}
    degree_v1 = degree_p + 1;
    degree_v2 = degree_v1;
    regularity_v = regularity_p;
    nsub_v = nsub_p;
    knots_v1 = kntrefine (knots, nsub_v-1, degree_v1, regularity_v);
    knots_v2 = knots_v1;

  case {'sg'}
    degree_v1 = degree_p + 1;
    degree_v2 = degree_v1;
    regularity_v = regularity_p+1;
    nsub_v = 2*nsub_p;
    knots_v1 = kntrefine (knots, nsub_v-1, degree_v1, regularity_v);
    knots_v2 = knots_v1;

  case {'ndl'}
    degree_v1 = degree_p+1;
    degree_v2 = degree_p+1;
    knots_v1{1}  = [knotsp{1}(1) knotsp{1} knotsp{1}(end)];
    knots_v1{2}  = sort ([knotsp{2}, unique(knotsp{2})]);
    knots_v2{1}  = sort ([knotsp{1}, unique(knotsp{1})]);
    knots_v2{2}  = [knotsp{2}(1) knotsp{2} knotsp{2}(end)];

  case {'rt'}
    degree_v1 = [degree_p(1)+1 degree_p(2)];
    degree_v2 = [degree_p(1) degree_p(2)+1];
    knots_v1{1}  = [knotsp{1}(1) knotsp{1} knotsp{1}(end)];
    knots_v1{2}  = knotsp{2};
    knots_v2{1}  = knotsp{1};
    knots_v2{2}  = [knotsp{2}(1) knotsp{2} knotsp{2}(end)];
  otherwise
    error ('sp_bspline_fluid_2d_phys: unknown element type')
end

% Construction of the space for the velocity
switch (lower (element_name))
  case {'th'}
    sp1 = sp_bspline_2d (knots_v1, degree_v1, msh);
    sp2 = sp_bspline_2d (knots_v2, degree_v2, msh);
    spv = sp_vector_2d (sp1, sp2, msh);

    PI = speye (spp.ndof);
  case {'sg'}
    sp1 = sp_bspline_2d (knots_v1, degree_v1, msh);
    sp2 = sp_bspline_2d (knots_v2, degree_v2, msh);
    spv = sp_vector_2d (sp1, sp2, msh);

    PI = speye (spp.ndof);
  case {'ndl'}
    sp1 = sp_bspline_2d (knots_v1, degree_v1, msh);
    sp2 = sp_bspline_2d (knots_v2, degree_v2, msh);
    spv = sp_vector_2d_piola_transform (sp1, sp2, msh);

    PI = speye (spp.ndof);
  case {'rt'}
    sp1 = sp_bspline_2d (knots_v1, degree_v1, msh);
    sp2 = sp_bspline_2d (knots_v2, degree_v2, msh);
    spv = sp_vector_2d_piola_transform (sp1, sp2, msh);

    PI = b2nst__ (spp, knotsp, degree_p, msh);
  otherwise
    error ('sp_bspline_fluid_2d_phys: unknown element type')
end

end

