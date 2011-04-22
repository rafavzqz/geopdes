% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011 Andrea Bressan, Rafael Vazquez
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

% This is fun_transform for 'TH' and 'SG' elements
function spv = do_spv_component_wise_2d__ (knotsv1, degreev1, knotsv2, degreev2, msh)
  spv1 = sp_bspline_2d_phys (knotsv1, degreev1, msh);
  spv2 = sp_bspline_2d_phys (knotsv2, degreev2, msh);
  spv  = sp_scalar_to_vector_2d (spv1, spv2, msh, 'divergence', true);
end