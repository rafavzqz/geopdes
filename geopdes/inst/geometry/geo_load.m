% GEO_LOAD: create a geometry structure from a file or a cell-array of functions.
%
% geometry = geo_load (in [, embedInR3])
%
% INPUT :
%
%   The input variable in may be either
%   - a structure representing a NURBS surface or volume, as in the NURBS toolbox
%   - a string variable with the name of the file to be read (see doc/geo_specs_v21.txt)
%   - a cell-array of function handles, to evaluate the function, the first order
%       derivatives, and eventually the second order derivatives
%   - a 4x4 matrix representing an affine transformation
%   - a gsTHBSpline object coming from G+smo
%   - a gsTensorBSpline object coming from G+smo
%   
%   embedInR3: boolean (default: false), true if the input geometry has to
%     be considered in R^3 in any case.
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
% Copyright (C) 2018 Ondine Chanon
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

function geometry = geo_load (in, embedInR3)

  if nargin < 2
    embedInR3 = false;
  end

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
    
%% geometry is given as a gsTHBSpline from G+smo
  elseif (isa (in, 'gsTHBSpline'))
    geometry.gismo = in;
    thsb = in.basis;
	geometry.map = @(ps) gismo_map (ps, in, 0); 
    geometry.map_der = @(ps) gismo_map (ps, in, 1); 
    geometry.map_der2 = @(ps) gismo_map (ps, in, 2);
    geometry.rdim = in.geoDim;
    geometry.order = zeros (1, in.parDim);
    geometry.regularity = zeros (thsb.maxLevel, in.parDim);
    
    for dir = 1:in.parDim
      orderDir = thsb.degree(dir) + 1;
      geometry.order(dir) = orderDir; 
      
      for lev = 1:thsb.maxLevel
        knotsLevDir = thsb.knots(lev, dir);
        geometry.knots{lev}{dir} = knotsLevDir;
        
        midknots = knotsLevDir(orderDir+1:end-orderDir);
        if isempty(midknots)
          geometry.regularity(lev, dir) = orderDir - 2;
        else
          geometry.regularity(lev, dir) = orderDir-1 - max (histc (midknots, ...
                                                       unique (midknots)));
        end
      end
    end
    
%% geometry is given as a gsTensorBSpline from G+smo
  elseif (isa (in, 'gsTensorBSpline'))
    geometry.gismo = in;
    geometry.map = @(ps) gismo_map (ps, in, 0); 
    geometry.map_der = @(ps) gismo_map (ps, in, 1); 
    geometry.map_der2 = @(ps) gismo_map (ps, in, 2);
    geometry.rdim = in.geoDim;
    geometry.order = zeros (1, in.parDim);
    geometry.regularity = zeros (1, in.parDim);
    
    thsb = in.basis;
    for dir = 1:in.parDim
      orderDir = thsb.degree(dir) + 1;
      geometry.order(dir) = orderDir; 
      
      knotsDir = thsb.knots(dir);
      geometry.knots{dir} = knotsDir;
      
      midknots = knotsDir(orderDir+1:end-orderDir);
	  if isempty(midknots)
		geometry.regularity(dir) = orderDir - 2;
	  else
		geometry.regularity(dir) = orderDir-1 - max (histc (midknots, ...
												     unique (midknots)));
	  end
    end
    
  else
    error ('geo_load: wrong input type');
  end
  
% In case the geometry is a NURBS structure, we use the NURBS toolbox
% This allows to use a rectangular parametric domain, instead of [0,1]^n.
  if (isfield (geometry, 'nurbs'))
    if ~embedInR3
      if (any (abs(geometry.nurbs.coefs(3,:)) > 1e-12))
        rdim = 3;
      elseif (any (abs(geometry.nurbs.coefs(2,:)) > 1e-12))
        rdim = 2;
      else
        rdim = 1;
      end
    else
      rdim = 3;
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
    
  elseif (isa (in, 'gsTHBSpline') || isa (in, 'gsTensorBSpline'))
    if in.parDim > 1
      for ibnd = 1:2*in.parDim
        geometry.boundary(ibnd).map     = @(PTS) boundary_map (geometry.map, ibnd, PTS);
        geometry.boundary(ibnd).map_der = @(PTS) boundary_map_der (geometry.map, geometry.map_der, ibnd, PTS);
      end
    else
      warning('Load of boundary map, map_der of G+smo geometry of ndim = 1 not implemented yet.')
    end
    
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
    ndim = size (pts,1) + 1;
    ind = setdiff (1:ndim, ind2);
    
    pts_aux(ind,:) = pts;
    if (mod (iside,2) == 0)
      pts_aux(ind2,:) = ones (1, size(pts,2));
    else
      pts_aux(ind2,:) = zeros (1, size(pts,2));
    end
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
      ndim = size (pts, 1) + 1;
      ind = setdiff (1:ndim, ind2);
      
      pts_aux(ind,:) = pts;
      if (mod (iside,2) == 0)
          pts_aux(ind2,:) = ones (1, size(pts,2));
      else
          pts_aux(ind2,:) = zeros (1, size(pts,2));
      end
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

% These two functions are used to compute mappings from G+smo geometries
function varargout = gismo_map (pts, in, der)
  if ( iscell (pts) ) 
      ndim = length (pts);
      npts = prod (cellfun (@length,pts));
      pts = cartesian_product_from_cell (pts);
  else
      [ndim, npts] = size (pts);
  end
  
  if (der == 0 || nargout > 1)
    F = in.eval(pts); % dimension rdim x npts
    varargout{1} = F;
  end
  if (der == 1 || nargout > 2)
    % g+smo dim: rdim x (ndim x npts)
    % geopdes dim: rdim x ndim x npts
    jac = reshape (in.jacobian(pts), [], ndim, npts); 
    if nargout == 1
        varargout{1} = jac;
    else
        varargout{2} = jac;
    end
  end
  if (der == 2)
    rdim = in.geoDim;
    % g+smo dim: rdim, (ndim x ndim) x npts
    % geopdes dim rdim x ndim x ndim x npts
    hess = zeros (rdim, ndim, ndim, npts);
    for dir = 1:rdim
      hess(dir,:,:,:) = reshape (in.hess(pts, dir), ndim, ndim, npts);
    end
    if nargout == 1
      varargout{1} = hess;
    elseif nargout == 3
      varargout{3} = hess;
    end
  end
end

function pts_aux = cartesian_product_from_cell (pts)
  % create cartesian product points from cell information
  s = cellfun (@length, pts);
  s_cell = cell (length(s), 1);
  for ii = 1:length(s)
      s_cell{ii} = 1:s(ii);
  end
  x = cell (1, numel (s_cell));
  [x{:}] = ndgrid (s_cell{:});
  pts_aux = [];
  for ii = 1:length(s)
      pts_aux = [pts_aux ; reshape( pts{ii}(x{ii}), 1, [] )];
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
