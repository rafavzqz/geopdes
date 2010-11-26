% SP_BSPLINE_TH_2D_PHYS: Construct the Taylor-Hood pair of B-Splines spaces on the physical domain in 2D.
%
%     [spv, spp] = sp_bspline_th_2d_phys (knotsp, degree, msh)
%
% INPUTS:
%
%     knotsp:       b-spline knots of the pressure space along each parametric direction
%     degree:       b-spline degree of the pressure space along each parametric direction
%     msh:          structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%
% OUTPUT:
%
%    spv: structure representing the discrete velocity function space
%    spp: structure representing the discrete pressure function space
%         each of the above has the following fields:
%
%
%   For more details, see:
%      A.Buffa, C.de Falco, G. Sangalli, 
%      IsoGeometric Analysis: Stable elements for the 2D Stokes equation
%      IJNMF, 2010
%
% Copyright (C) 2009, 2010 Carlo de Falco
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

function [spv, spp] = sp_bspline_th_2d_phys (knotsp, degree, msh)

for idim = 1:2
  knotsv{idim}  = sort ([knotsp{idim}, unique(knotsp{idim})]);
end
spp  = sp_bspline_2d_phys (knotsp, degree, msh, 'gradient', false);

spv = do_spv (knotsv, degree, msh);
spv.spfun  = @(MSH) do_spv (knotsv, degree, MSH);
end

function spv = do_spv (knotsv, degree, msh)
  spv  = sp_bspline_2d_phys (knotsv, degree+1, msh);
  spv  = sp_scalar_to_vector_2d (spv, spv, msh, 'divergence', true);
end
