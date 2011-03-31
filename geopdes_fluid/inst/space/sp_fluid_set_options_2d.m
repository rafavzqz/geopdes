% SP_FLUID_SET_OPTIONS_2D: Construct the knot vectors and set some options for the construction of different pair of spaces for fluid problems.
%
% 
% [knotsp, knotsv1, degreev1, knotsv2, degreev2, fun_trans, press_proj, der2] = 
%  sp_fluid_set_options (elem_name, knots, nsub_p, degree_p, regularity_p)
%
% INPUTS:
%
%   elem_name:    the type of element. Available choices are 'TH' (Taylor-Hood),
%                  'NDL' (Nedelec, 2nd family), 'RT' (Raviart-Thomas) and 
%                  'SG' (subgrid). For more details see the references below
%   knots:        knot vector to be refined (usually, the one of the geometry)
%   nsub_p:       number of subdivisions O QUALCOSA DEL GENERE. DEVO SISTEMARE
%   degree_p:     degree of the pressure space to be computed later
%   regularity_p: regularity of the spline pressure space to be computed later
%
% OUTPUT:
%
%   knotsp:     knot vector of the pressure space along each parametric direction
%   knotsv1:    knot vector of the space for the first parametric component of 
%                the velocity along each parametric direction
%   degreev1:   degree of the space for the first parametric component of
%                the velocity along each parametric direction
%   knotsv2:    knot vector of the space for the second parametric component of 
%                the velocity along each parametric direction
%   degreev2:   degree of the space for the second parametric component of
%                the velocity along each parametric direction
%   fun_trans:  transformation of the velocity from the parametric to the
%                physical domain. Component-wise for TH-SG, Piola for NDL-RT.
%   der2:       option to know whether the second derivative of the
%                parameterization must be computed or not
%   press_proj: function to compute a projector for the pressure, needed
%                for RT spaces
%
%   For more details, see:
%      A.Buffa, C.de Falco, G. Sangalli, 
%      IsoGeometric Analysis: Stable elements for the 2D Stokes equation
%      IJNMF, 2010
%
% Copyright (C) 2011 Andrea Bressan, Carlo de Falco, Rafael Vazquez
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

function [knotsp, knotsv1, degreev1, knotsv2, degreev2, der2] = ...
    sp_fluid_set_options_2d (elem_name, knots, nsub_p, degree_p, regularity_p)

  knotsp = kntrefine (knots, nsub_p-1, degree_p, regularity_p);

  switch (lower (elem_name))
   case 'th'
    der2 = false;

    degreev1 = degree_p + 1;
    degreev2 = degreev1;
    regularity_v = regularity_p;
    nsub_v = nsub_p;
    knotsv1 = kntrefine (knots, nsub_v-1, degreev1, regularity_v);
    knotsv2 = knotsv1;

%    fun_transform = @(knotsv1, degreev1, knotsv2, degreev2, msh) ...
%        do_spv_component_wise__ (knotsv1, degreev1, knotsv2, degreev2, msh);
%    press_proj = @(spp, knotsp, degreep, msh) speye (spp.ndof);
   case 'sg'
    der2 = false;

    degreev1 = degree_p + 1;
    degreev2 = degreev1;
    regularity_v = regularity_p+1;
    nsub_v = 2*nsub_p;
    knotsv1 = kntrefine (knots, nsub_v-1, degreev1, regularity_v);
    knotsv2 = knotsv1;

%    fun_transform = @(knotsv1, degreev1, knotsv2, degreev2, msh) ...
%        do_spv_component_wise__ (knotsv1, degreev1, knotsv2, degreev2, msh);
%    press_proj = @(spp, knotsp, degreep, msh) speye (spp.ndof);
   case 'ndl'
    der2 = true;

    knotsv1{1}  = [knotsp{1}(1) knotsp{1} knotsp{1}(end)];
    knotsv1{2}  = sort ([knotsp{2}, unique(knotsp{2})]);
    knotsv2{1}  = sort ([knotsp{1}, unique(knotsp{1})]);
    knotsv2{2}  = [knotsp{2}(1) knotsp{2} knotsp{2}(end)];
    degreev1 = degree_p+1;
    degreev2 = degree_p+1;

%    fun_transform = @(knotsv1, degreev1, knotsv2, degreev2, msh) ...
%        do_spv_piola_transform__ (knotsv1, degreev1, knotsv2, degreev2, msh);
%    press_proj = @(spp, knotsp, degreep, msh) speye (spp.ndof);
   case 'rt'
    der2 = true;

    knotsv1{1}  = [knotsp{1}(1) knotsp{1} knotsp{1}(end)];
    knotsv1{2}  = knotsp{2};
    knotsv2{1}  = knotsp{1};
    knotsv2{2}  = [knotsp{2}(1) knotsp{2} knotsp{2}(end)];
    degreev1 = [degree_p(1)+1 degree_p(2)];
    degreev2 = [degree_p(1) degree_p(2)+1];

%    fun_transform = @(knotsv1, degreev1, knotsv2, degreev2, msh) ...
%        do_spv_piola_transform__ (knotsv1, degreev1, knotsv2, degreev2, msh);
%    press_proj = @(spp, knotsp, degreep, msh) b2nst_odd__ (spp, knotsp, degreep, msh);
   otherwise
    error('space_set: Unknown element type')
  end

end
