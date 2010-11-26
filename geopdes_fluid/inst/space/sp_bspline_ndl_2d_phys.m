% SP_BSPLINE_NDL_2D_PHYS: Construct the Nedelec's second family pair of B-Splines spaces on the physical domain in 2D.
%
%     [spv, spp] = sp_bspline_ndl_2d_phys (knotsp, degree, msh)
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

function [spv, spp] = sp_bspline_ndl_2d_phys (knotsp, degree, msh)

knotsvx{1}  = [knotsp{1}(1) knotsp{1} knotsp{1}(end)];
knotsvx{2}  = sort ([knotsp{2}, unique(knotsp{2})]);
knotsvy{1}  = sort ([knotsp{1}, unique(knotsp{1})]);
knotsvy{2}  = [knotsp{2}(1) knotsp{2} knotsp{2}(end)];
dvx = degree+1;
dvy = degree+1;


spp = sp_bspline_2d_phys (knotsp, degree, msh, 'gradient', false);

spv = do_spv__ (knotsvx, dvx, knotsvy, dvy, msh);
spv.spfun = @(MSH) do_spv__ (knotsvx, dvx, knotsvy, dvy, MSH);

end

function spv = do_spv__ (knotsvx, dvx, knotsvy, dvy, msh);
spvx  = sp_bspline_2d_param (knotsvx, dvx, msh);
spvy  = sp_bspline_2d_param (knotsvy, dvy, msh);
spv   = sp_scalar_to_vector_2d (spvx, spvy, msh, 'divergence', true);
spv   = sp_piola_transform_2d  (spv, msh);
end