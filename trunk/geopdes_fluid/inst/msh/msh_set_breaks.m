% MSH_SET_BREAKS: Construct the breaks for the quadrature rule, depending on the choice of the discrete space.
%
% [breaks, der2] = msh_set_breaks (elem_name, knots, nsub)
%
% INPUTS:
%
%   elem_name:    the type of element. Available choices are 'TH' (Taylor-Hood),
%                  'NDL' (Nedelec, 2nd family), 'RT' (Raviart-Thomas) and 
%                  'SG' (subgrid). For more details see the references below
%   knots:        knot vector of the coarse geometry
%   nsub:         number of subdivisions for refinement of the pressure space.
%
% OUTPUT:
%
%   breaks:     interval breaks for the quadrature rule. 
%   der2:       option to know whether the second derivative of the
%                parameterization must be computed or not
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
% Copyright (C) 2011 Andrea Bressan, Carlo de Falco, Rafael Vazquez
% Copyright (C) 2014 Elena Bulgarello, Carlo de Falco, Sara Frizziero
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

function [breaks, der2] = msh_set_breaks (elem_name, knots, nsub)

  ndim = numel (knots);
  degree1 = ones (ndim, 1);
  reg0 = zeros (ndim, 1);

  switch (lower (elem_name))
   case {'th'}
    der2 = false;
    breaks = kntrefine (knots, nsub-1, degree1, reg0);

   case {'sg'}
    der2 = false;
    breaks = kntrefine (knots, 2*nsub-1, degree1, reg0);

   case {'ndl'}
    if (numel (knots) == 2)
      der2 = true;
      breaks = kntrefine (knots, nsub-1, degree1, reg0);
    else
      error ('NDL elements have not been implemented in 3D yet')
    end
   case {'rt'}
     der2 = true;
     breaks = kntrefine (knots, nsub-1, degree1, reg0);
   otherwise
    error('msh_set_breaks: Unknown element type')
  end

end
