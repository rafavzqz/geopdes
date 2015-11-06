% GEO_LOAD: create a geometry structure from a file or a cell-array of functions.
%
% geometry = geo_load (input)
%
% INPUT :
%
%   the input variable may be either
%   - a string variable with the name of the file to be read
%   - a cell-array of function handles
%   - a 4x4 matrix representing an affine transformation
%   - a structure representing a nurbs surface or volume
%
% OUTPUT:
%
%   geometry: a structure that contains, at least, the following fields
%             map:      a function handle to evaluate the parameterization
%             map_der:  a function handle to evaluate the derivatives of the parameterization
%             map_der2: a function handle to evaluate the second derivatives of the parameterization
%   The structure may contain further information. See the documentation.
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2013 Rafael Vazquez
% Copyright (C) 2014 Elena Bulgarello, Carlo de Falco, Sara Frizziero
% Copyright (C) 2015 Rafael Vazquez
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

function geometry = geo_load (in)

  if (isstruct (in) && isfield (in, 'form') && strcmpi (in.form, 'B-NURBS')) 
%% geometry is given as a NURBS struct
    geometry.nurbs   = in;

  elseif (ischar (in))
%% load geometry from a file
    if (strcmpi (in(end-3:end), '.mat'))     
      tmp  = load (in);
      geometry.nurbs   = tmp.geo;
    elseif (strcmpi (in(end-3:end), '.txt'))
      geometry = geo_read_nurbs (in);
%% load geometry from an xml file
    elseif (strcmpi (in(end-3:end), '.xml'))
      geometry.nurbs = xml_import (in);
    else
      error ('geo_load: unknown file extension');
    end

  elseif (iscell (in))
%% geometry is given as a cell array of functions
    geometry.map     =  in{1};
    geometry.map_der =  in{2};
    if (length (in) > 2)
      geometry.map_der2 =  in{3};
    end
    
  elseif (isnumeric (in) && all (size (in) == [4, 4]))  
%% geometry is given as a 4x4 matrix representing an affine transformation
    geometry.map = @(ps) affine_map  (ps, in);
    geometry.map_der = @(ps) affine_map_der  (ps, in);
    geometry.map_der2 = @affine_map_der2;
  else
    error ('geo_load: wrong input type');
  end
  
% In case the geometry is a NURBS structure, we use the NURBS toolbox
% This allows to use a rectangular parametric domain, instead of [0,1]^n.
  if (isfield (geometry, 'nurbs'))
    if (any (abs(geometry.nurbs.coefs(3,:)) > 1e-12))
      rdim = 3;
    elseif (any (abs(geometry.nurbs.coefs(2,:)) > 1e-12))
      rdim = 2;
    else
      rdim = 1;
    end
    geometry.rdim = rdim;

    warning ('off', 'nrbderiv:SecondDerivative')
    [deriv, deriv2] = nrbderiv (geometry.nurbs);
    geometry.dnurbs = deriv;
    geometry.dnurbs2 = deriv2;

    geometry.map      =  @(PTS) geo_nurbs (geometry.nurbs, deriv, deriv2, PTS, 0, rdim);
    geometry.map_der  =  @(PTS) geo_nurbs (geometry.nurbs, deriv, deriv2, PTS, 1, rdim);
    geometry.map_der2 =  @(PTS) geo_nurbs (geometry.nurbs, deriv, deriv2, PTS, 2, rdim);

    if (numel (geometry.nurbs.order) > 1)
      bnd = nrbextract (geometry.nurbs);
      for ibnd = 1:numel (bnd)
        [deriv, deriv2] = nrbderiv (bnd(ibnd));
        geometry.boundary(ibnd).nurbs    = bnd(ibnd);
        geometry.boundary(ibnd).dnurbs   = deriv;
        geometry.boundary(ibnd).dnurbs2  = deriv2;
        geometry.boundary(ibnd).rdim     = rdim;
        geometry.boundary(ibnd).map      = @(PTS) geo_nurbs (bnd(ibnd), deriv, deriv2, PTS, 0, rdim);
        geometry.boundary(ibnd).map_der  = @(PTS) geo_nurbs (bnd(ibnd), deriv, deriv2, PTS, 1, rdim);
        geometry.boundary(ibnd).map_der2 = @(PTS) geo_nurbs (bnd(ibnd), deriv, deriv2, PTS, 2, rdim);
      end
    end
    warning ('on', 'nrbderiv:SecondDerivative')
  else
    for ibnd = 1:6 % This loop should be until 2*ndim, but ndim is not known
      geometry.boundary(ibnd).map     = @(PTS) boundary_map (geometry.map, ibnd, PTS);
      geometry.boundary(ibnd).map_der = @(PTS) boundary_map_der (geometry.map, geometry.map_der, ibnd, PTS);
    end
  end

end

function mps = affine_map  (ps, in)
  if (iscell (ps))
    ndim = numel (ps);
    npts = cellfun (@numel, ps);
    nps = prod (npts);
    if (ndim == 2)
      u = reshape (repmat (ps{1}(:), 1, npts(2)), 1, []);
      v = reshape (repmat (ps{2}(:)', npts(1), 1), 1, []);
      ps = [u; v];
    elseif (ndim == 3)
      u = reshape (ps{1}, npts(1), 1, 1);
      u = reshape (repmat (u, [1, npts(2), npts(3)]), 1, []);
      v = reshape (ps{2}, 1, npts(2), 1);
      v = reshape (repmat (v, [npts(1), 1, npts(3)]), 1, []);
      w = reshape (ps{3}, 1, 1, npts(3));
      w = reshape (repmat (w, [npts(1), npts(2), 1]), 1, []);
      ps = [u; v; w];
    end
  else
    ndim = size (ps, 1);
    nps  = size (ps, 2);
  end
  mps  = in([1:ndim, 4], [1:ndim, 4]) * [ps; ones(1, nps)];
  mps  = mps (1:ndim, :);
end

function mps = affine_map_der  (ps, in)
  if (iscell (ps))
    ndim = numel (ps);
    npts = cellfun (@numel, ps);
    nps = prod (npts);
    if (ndim == 2)
      u = reshape (repmat (ps{1}(:), 1, npts(2)), 1, []);
      v = reshape (repmat (ps{2}(:)', npts(1), 1), 1, []);
      ps = [u; v];
    elseif (ndim == 3)
      u = reshape (ps{1}, npts(1), 1, 1);
      u = reshape (repmat (u, [1, npts(2), npts(3)]), 1, []);
      v = reshape (ps{2}, 1, npts(2), 1);
      v = reshape (repmat (v, [npts(1), 1, npts(3)]), 1, []);
      w = reshape (ps{3}, 1, 1, npts(3));
      w = reshape (repmat (w, [npts(1), npts(2), 1]), 1, []);
      ps = [u; v; w];
    end
  else
    ndim = size (ps, 1);
    nps  = size (ps, 2);
  end
  mps  = repmat (in(1:ndim, 1:ndim), [1, 1, nps]);
end

function mps = affine_map_der2  (ps)
  if (iscell (ps))
    ndim = numel (ps);
    nps = prod (cellfun (@numel, ps));
  else
    ndim = size (ps, 1);
    nps  = size (ps, 2);
  end
  mps  = zeros (ndim, ndim, ndim, nps);
end


% These two functions are to compute boundary entities from the global ones
function F = boundary_map (map, iside, pts)

%%    ind  = [2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind  = [2 2 1 1] in 2D;
%%    ind2 = [1 1 2 2 3 3] in 3D,                  %ind2 = [1 1 2 2] in 2D
  ind2 = ceil (iside/2);
  
  if (iscell (pts))
    ndim = numel (pts) + 1;
    ind = setdiff (1:ndim, ind2);

    pts_aux(ind) = pts;
    if (mod (iside, 2) == 1)
      pts_aux{ind2} = 0;
    else
      pts_aux{ind2} = 1;
    end
  else
    error ('For the boundary, a cell array should be passed as the argument')
  end

  F = map (pts_aux);

end

function varargout = boundary_map_der (map, map_der, iside, pts)

%%    ind  = [2 3; 2 3; 1 3; 1 3; 1 2; 1 2] in 3D, %ind  = [2 2 1 1] in 2D;
%%    ind2 = [1 1 2 2 3 3] in 3D,                  %ind2 = [1 1 2 2] in 2D
  ind2 = ceil (iside/2);
  
  if (iscell (pts))
    ndim = numel (pts) + 1;
    ind = setdiff (1:ndim, ind2);

    pts_aux(ind) = pts;
    if (mod (iside, 2) == 1)
      pts_aux{ind2} = 0;
    else
      pts_aux{ind2} = 1;
    end
  else
    error ('For the boundary, a cell array should be passed as the argument')
  end

  DF = map_der (pts_aux);
  DF = DF(:,ind,:);
  if (nargout == 1)
    varargout{1} = DF;
  elseif (nargout == 2)
    F = map (pts_aux);
    varargout{1} = F;
    varargout{2} = DF;
  end

end


%!shared g1,g2,x
%!test
%! g1 = geo_load ('geo_ring.txt');
%! g2 = geo_load ('ring.mat');
%! for [v, k] = g1.nurbs
%!  assert (isequal (v, g2.nurbs.(k)))
%! endfor
%! x = rand (2, 10);
%!test
%! assert (g1.map (x), g2.map (x))
%!test
%! assert (g1.map_der (x), g2.map_der (x))
