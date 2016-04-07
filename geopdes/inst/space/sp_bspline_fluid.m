% SP_BSPLINE_FLUID: Construct different pair of B-Splines spaces on the physical domain for fluid problems.
%
%   [spv, spp] = sp_bspline_fluid (elem_name, knots, nsub, ...
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
%                or points for visualization (see msh_cartesian).
%
% OUTPUT:
%
%   spv: object representing the discrete velocity function space (see sp_vector)
%   spp: object representing the discrete pressure function space (see sp_scalar)
%
%   For more details, see:
%      A. Buffa, C. de Falco, G. Sangalli, 
%      IsoGeometric Analysis: Stable elements for the 2D Stokes equation
%      Internat. J. Numer. Methods Fluids, 2011
%
%      A. Bressan, G. Sangalli,
%      Isogeometric discretizations of the Stokes problem: stability
%       analysis by the macroelement technique
%      IMA J. Numer. Anal., 2013.
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, 2015 Rafael Vazquez
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

function [spv, spp] = sp_bspline_fluid (element_name, ...
                   knots, nsub_p, degree_p, regularity_p, msh)

% Construction of the knot vectors and the discrete space for the velocity
switch (lower (element_name))
  case {'th'}
    knotsp = kntrefine (knots, nsub_p-1, degree_p, regularity_p);
    spp = sp_bspline (knotsp, degree_p, msh);
  
    degree_v = degree_p + 1;
    regularity_v = regularity_p;
    nsub_v = nsub_p;
    knots_v = kntrefine (knots, nsub_v-1, degree_v, regularity_v);
    scalar_space = sp_bspline (knots_v, degree_v, msh);
    for idim = 1:msh.ndim
      scalar_spaces{idim} = scalar_space;
    end
    spv = sp_vector (scalar_spaces, msh);
    clear scalar_spaces scalar_space

  case {'sg'}
    knotsp = kntrefine (knots, nsub_p-1, degree_p, regularity_p);
    spp = sp_bspline (knotsp, degree_p, msh);
    
    degree_v = degree_p + 1;
    regularity_v = regularity_p+1;
    nsub_v = 2*nsub_p;
    knots_v = kntrefine (knots, nsub_v-1, degree_v, regularity_v);
    scalar_space = sp_bspline (knots_v, degree_v, msh);
    for idim = 1:msh.ndim
      scalar_spaces{idim} = scalar_space;
    end
    spv = sp_vector (scalar_spaces, msh);
    clear scalar_spaces scalar_space

  case {'ndl'}
% In this case the regularity is assigned first in the velocity space
    degree_h1 = degree_p + 1;
    regularity_h1 = regularity_p + 1;
    knots_h1 = kntrefine (knots, nsub_p-1, degree_h1, regularity_h1);
    [knots_hdiv, degree_hdiv] = knt_derham (knots_h1, degree_h1, 'Hdiv');

    degree_v = degree_h1;
    for idim = 1:msh.ndim
      knots_v{idim} = knots_hdiv{idim}{idim};
      for jdim = setdiff (1:msh.ndim, idim)
        knots_v{jdim} = sort ([knots_hdiv{idim}{jdim}, unique(knots_hdiv{idim}{jdim})]);
      end
      scalar_spaces{idim} = sp_bspline (knots_v, degree_v, msh);
    end
    spv = sp_vector (scalar_spaces, msh, 'div-preserving');
    clear scalar_spaces
    
    [knotsp, degp] = knt_derham (knots_h1, degree_h1, 'L2');
    spp = sp_bspline (knotsp, degp, msh, 'integral-preserving');

  case {'rt'}
% In this case the regularity is assigned first in the velocity space
    degree_h1 = degree_p + 1;
    regularity_h1 = regularity_p + 1;
    knots_h1 = kntrefine (knots, nsub_p-1, degree_h1, regularity_h1);

    [knots_v, degree_v] = knt_derham (knots_h1, degree_h1, 'Hdiv');
    for idim = 1:msh.ndim
      scalar_spaces{idim} = sp_bspline (knots_v{idim}, degree_v{idim}, msh);
    end
    spv = sp_vector (scalar_spaces, msh, 'div-preserving');
    clear scalar_spaces

    [knotsp, degp] = knt_derham (knots_h1, degree_h1, 'L2');
    spp = sp_bspline (knotsp, degp, msh, 'integral-preserving');

%     if (nargout == 3)
%       PI = b2nst__ (spp, knotsp, degree_p, msh);
%     end
  otherwise
    error ('sp_bspline_fluid: unknown element type')
end

end

