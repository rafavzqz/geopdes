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
%             map:     a function handle to evaluate the parameterization
%             map_der: a function handle to evaluate the derivatives of the parameterization
%   The structure may contain further information. See the documentation.
%
% Copyright (C) 2010 Carlo de Falco
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
%
% Author: Carlo de Falco <cdf@users.sourceforge.net>
% Created: 2010-05-06

function geometry = geo_load (in)

  if (isstruct (in) && isfield (in, 'form') && strcmpi (in.form, 'B-NURBS')) 
%% geometry is given as a NURBS struct
    geometry.nurbs   = in;
    if (numel (geometry.nurbs.knots) == 2)
      geometry.map      =  @(PTS) geo_2d_nurbs (geometry.nurbs, PTS, 0);
      geometry.map_der  =  @(PTS) geo_2d_nurbs (geometry.nurbs, PTS, 1);
      geometry.map_der2 =  @(PTS) geo_2d_nurbs (geometry.nurbs, PTS, 2);
    elseif (numel (geometry.nurbs.knots) == 3)
      geometry.map     =  @(PTS) geo_3d_nurbs (geometry.nurbs, PTS, 0);
      geometry.map_der =  @(PTS) geo_3d_nurbs (geometry.nurbs, PTS, 1);
    end

  elseif (ischar (in))
%% load geometry from a file
    if (strcmpi (in(end-3:end), '.mat'))     
      tmp  = load (in);
      geometry.nurbs   = tmp.geo;
      if (numel (geometry.nurbs.knots) == 2)
        geometry.map      =  @(PTS) geo_2d_nurbs (geometry.nurbs, PTS, 0);
        geometry.map_der  =  @(PTS) geo_2d_nurbs (geometry.nurbs, PTS, 1);
        geometry.map_der2 =  @(PTS) geo_2d_nurbs (geometry.nurbs, PTS, 2);
    elseif (numel (geometry.nurbs.knots) == 3)
        geometry.map     =  @(PTS) geo_3d_nurbs (geometry.nurbs, PTS, 0);
        geometry.map_der =  @(PTS) geo_3d_nurbs (geometry.nurbs, PTS, 1);
      end
    elseif (strcmpi (in(end-3:end), '.txt'))
      geometry = geo_read_nurbs (in);
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
%% geometry is given as a 4x4 matrix representig an affine transformation
    geometry.map = @(ps) affine_map  (ps, in);
    geometry.map_der = @(ps) affine_map_der  (ps, in);
    geometry.map_der2 = @affine_map_der2;
  else
    error ('geo_load: wrong input type');
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
