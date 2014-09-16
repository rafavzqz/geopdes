% SP_BSPLINE_FLUID_3D: Construct different pair of B-Splines spaces on the physical domain for fluid problems.
%
%   [spv, spp] = sp_bspline_fluid_3d (elem_name, knots, nsub, ...
%                                              degreep, regularity, msh)
%
% INPUTS:
%
%   elem_name: the name of the element. Right now 'TH' (Taylor-Hood), 
%              and 'SG' (SubGrid) are supported.
%   knots:      knot vector of the coarse geometry.
%   nsub:       number of subdivisions of each interval.
%   degreep:    degree of the pressure space along each parametric direction
%   regularity: continuity of the pressure space along each parametric direction
%   msh:        msh object containing (in the field msh.qn) the points 
%                along each parametric direction in the parametric 
%                domain at which to evaluate, i.e. quadrature points 
%                or points for visualization (see msh_3d).
%
% OUTPUT:
%
%   spv: object representing the discrete velocity function space (see sp_vector_3d)
%   spp: object representing the discrete pressure function space (see sp_bspline_3d)
%
%   For more details, see:
%      A. Buffa, C. de Falco, G. Sangalli, 
%      IsoGeometric Analysis: Stable elements for the 2D Stokes equation
%      Internat. J. Numer. Methods Fluids, 2011
%
%      A. Bressan, G. Sangalli,
%      Isogeometric discretizations of the Stokes problem: stability
%       analysis by the macroelement technique
%      Tech. Report, IMATI-CNR, 2011.
%
% Copyright (C) 2009, 2010, 2011, 2014 Carlo de Falco
% Copyright (C) 2011 Andrea Bressan, Rafael Vazquez
% Copyright (C) 2014 Elena Bulgarello, Sara Frizziero
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

function [spv, spp] = sp_bspline_fluid_3d (element_name, ...
                                           knots, nsub_p, degree_p, regularity_p, msh)

  %% The pressure space is the same in all the cases
  knotsp = kntrefine (knots, nsub_p-1, degree_p, regularity_p);
  spp = sp_bspline_3d (knotsp, degree_p, msh);

  switch (lower (element_name))
    case {'th'}
      degree_v = degree_p + 1;
      regularity_v = regularity_p;
      nsub_v = nsub_p;
      knots_v = kntrefine (knots, nsub_v-1, degree_v, regularity_v);
      sp_scalar = sp_bspline_3d (knots_v, degree_v, msh);
      spv = sp_vector_3d (sp_scalar, sp_scalar, sp_scalar, msh);

    case {'sg'}
      degree_v = degree_p + 1;
      regularity_v = regularity_p+1;
      nsub_v = 2*nsub_p;
      knots_v = kntrefine (knots, nsub_v-1, degree_v, regularity_v);
      sp_scalar = sp_bspline_3d (knots_v, degree_v, msh);
      spv = sp_vector_3d (sp_scalar, sp_scalar, sp_scalar, msh);

    case {'rt'}
      degree_v1 = [degree_p(1)+1 degree_p(2) degree_p(3)];
      degree_v2 = [degree_p(1) degree_p(2)+1 degree_p(3)];
      degree_v3 = [degree_p(1) degree_p(2) degree_p(3)+1];
      knots_v1{1}  = [knotsp{1}(1) knotsp{1} knotsp{1}(end)];
      knots_v1{2}  = knotsp{2};
      knots_v1{3}  = knotsp{3};
      knots_v2{1}  = knotsp{1};
      knots_v2{2}  = [knotsp{2}(1) knotsp{2} knotsp{2}(end)];
      knots_v2{3}  = knotsp{3};
      knots_v3{1}  = knotsp{1};
      knots_v3{2}  = knotsp{2};
      knots_v3{3}  = [knotsp{3}(1) knotsp{3} knotsp{3}(end)];

      sp1 = sp_bspline_3d (knots_v1, degree_v1, msh);
      sp2 = sp_bspline_3d (knots_v2, degree_v2, msh);
      sp3 = sp_bspline_3d (knots_v3, degree_v3, msh);
      spv = sp_vector_3d_piola_transform (sp1, sp2, sp3, msh);
    case {'ndl'}
      error ('NDL elements have not been implemented in 3D yet')
    otherwise
      error ('sp_bspline_fluid_3d: unknown element type')
  end

end
