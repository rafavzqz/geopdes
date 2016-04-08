% MP_GEO_LOAD: create a multipatch geometry structure and the interfaces structure from a file.
%
% [geometry, boundaries, interfaces] = mp_geo_load (input)
%
% INPUT :
%
%   input: a string variable with the name of the file to be read
%
% OUTPUT:
%
%   geometry:   an array of structures that contain the following fields
%               map:     a function handle to evaluate the parametrization
%               map_der: a function handle to evaluate the derivatives of the parametrization
%               nurbs:   a structure compatible with the NURBS toolbox
%
%   boundaries: an array of structures that contain the following fields
%               nsides:  number of faces conforming each boundary
%               patches: patches to which these faces belong
%               faces:   local numbering of the faces in their patch
%
%   interfaces: an array of structures that contain the following fields
%
%               patch1: the index of the first coinciding patch, 
%               side1:  the local index in patch1 of the interface boundary vertex (1D), edge (2D) or face (3D)
%               patch2: the index of the second coinciding patch
%               side2:  the local index in patch2 of the interface boundary vertex (1D), edge (2D) or face (3D)
%      2D Case
%               ornt:          a flag telling if the direction on patch1 matches the direction on patch2
%      3D Case
%               flag:          a flag telling if the u- and v- parameter directions on the first face coincide 
%                                   with the u- and v- directions on the second face
%               ornt1:         a flag telling if the u- direction on the first face matches the direction on patch2
%               ornt2:         a flag telling if the v- direction on the first face matches the direction on patch2
%
%   subdomains: an array of structures that contains the following fields
%               name:    the name of the subdomain
%               patches: indices of the patches belonging to the subdomain
%
%   boundary_interfaces: similar to interfaces, but for the (n-1) dimensional boundary
%
% This format is based on the paper
% T. Dokken, E. Quak, V. Skytt, 
%  Requirements from Isogeometric Analysis for Changes in Product Design Ontologies. 
%  Proceedings of the Focus K3D Conference on Semantic 3D Media and Content, Sophia Antipolis (France), 2010.
%
% Copyright (C) 2010, 2011 Carlo de Falco, Rafael Vazquez
% Copyright (C) 2013, 2015 Rafael Vazquez
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

function [geometry, boundaries, interfaces, subdomains, boundary_interfaces] = mp_geo_load (in)
  
  if (ischar (in))
    if (strcmpi (in(end-3:end), '.mat'))
%% load geometry from a file
      tmp = load(in);
      if (isfield (tmp, 'geometry') && isfield (tmp, 'boundaries') && isfield (tmp,'interfaces'))
        geometry = tmp.geometry;
        boundaries = tmp.boundaries;
        interfaces = tmp.interfaces;
        dim = numel (geometry(1).nurbs.knots);
        if (isfield (tmp, 'subdomains'))
          subdomains = tmp.subdomains;
        end
      elseif (isfield (tmp,'geo'))
        geometry.nurbs = tmp.geo;
        dim = numel (geometry.nurbs.knots);
        interfaces = [];
        if (isfield (tmp, 'boundaries'))
          boundaries = tmp.boundaries;
        else
          for iface = 1:2*dim
            boundaries(iface).name = ['BOUNDARY ' num2str(iface)];
            boundaries(iface).nsides = 1;
            boundaries(iface).patches = 1;
            boundaries(iface).faces = iface;
          end
        end
        if (isfield (tmp, 'subdomains'))
          subdomains = tmp.subdomains;
        else
          subdomains(1).name = 'SUBDOMAIN 1';
          subdomains(1).patches = 1;
        end
      end
      npatch = numel (geometry);

      for iptc = 1:npatch
        geometry(iptc).map      =  ...
          @(PTS) geo_nurbs (geometry(iptc).nurbs, PTS, 0);
        geometry(iptc).map_der  =  ...
          @(PTS) geo_nurbs (geometry(iptc).nurbs, PTS, 1);
        geometry(iptc).map_der2 =  ...
          @(PTS) geo_nurbs (geometry(iptc).nurbs, PTS, 2);
      end

    elseif (strcmpi (in(end-3:end), '.txt'))
%% load geometry from a txt file
      [geometry, boundaries, interfaces, subdomains] = mp_geo_read_nurbs (in);
    elseif (strcmpi (in(end-3:end), '.xml'))
%% load geometry from an xml file
      [nurbs, boundaries, interfaces] = xml_import (in);
      dim = numel (nurbs(1).knots);
      for iptc = 1:numel(nurbs)
        geometry(iptc).nurbs = nurbs(iptc);
        geometry(iptc).map      =  ...
          @(PTS) geo_nurbs (geometry(iptc).nurbs, PTS, 0);
        geometry(iptc).map_der  =  ...
          @(PTS) geo_nurbs (geometry(iptc).nurbs, PTS, 1);
        geometry(iptc).map_der2 =  ...
          @(PTS) geo_nurbs (geometry(iptc).nurbs, PTS, 2);
      end
    else
      error ('mp_geo_load: unknown file extension');
    end

  elseif (isstruct (in) && isfield (in, 'form') && strcmpi (in(1).form, 'B-NURBS'))
    for iptc = 1:numel(in)
      geometry(iptc) = geo_load (in(iptc));
    end
    warning ('Automatically generating the interface and boundary information with nrbmultipatch')
    [interfaces, boundaries] = nrbmultipatch (in);
    subdomains(1).name = 'SUBDOMAIN 1';
    subdomains(1).patches = 1:numel(in);
  else
    error ('mp_geo_load: wrong input type');
  end

  
  if (isfield (geometry, 'nurbs'))
    rdim = 0;
    for iptc = 1:numel (geometry)
      if (any (abs(geometry(iptc).nurbs.coefs(3,:)) > 1e-12))
        rdim = 3;
      elseif (any (abs(geometry(iptc).nurbs.coefs(2,:)) > 1e-12))
        rdim = max (rdim, 2);
      else
        rdim = max (rdim, 1);
      end
    end
    
    for iptc = 1:numel (geometry)
      geometry(iptc).rdim = rdim;
      
      [deriv, deriv2] = nrbderiv (geometry(iptc).nurbs);
      geometry(iptc).dnurbs = deriv;
      geometry(iptc).dnurbs2 = deriv2;

      geometry(iptc).map      =  @(PTS) geo_nurbs (geometry(iptc).nurbs, deriv, deriv2, PTS, 0, rdim);
      geometry(iptc).map_der  =  @(PTS) geo_nurbs (geometry(iptc).nurbs, deriv, deriv2, PTS, 1, rdim);
      geometry(iptc).map_der2 =  @(PTS) geo_nurbs (geometry(iptc).nurbs, deriv, deriv2, PTS, 2, rdim);

      if (numel (geometry(1).nurbs.order) > 1)
        bnd = nrbextract (geometry(iptc).nurbs);
        for ibnd = 1:numel (bnd)
          [deriv, deriv2] = nrbderiv (bnd(ibnd));
          geometry(iptc).boundary(ibnd).nurbs    = bnd(ibnd);
          geometry(iptc).boundary(ibnd).dnurbs   = deriv;
          geometry(iptc).boundary(ibnd).dnurbs2  = deriv2;
          geometry(iptc).boundary(ibnd).rdim     = rdim;
          geometry(iptc).boundary(ibnd).map      = @(PTS) geo_nurbs (bnd(ibnd), deriv, deriv2, PTS, 0, rdim);
          geometry(iptc).boundary(ibnd).map_der  = @(PTS) geo_nurbs (bnd(ibnd), deriv, deriv2, PTS, 1, rdim);
          geometry(iptc).boundary(ibnd).map_der2 = @(PTS) geo_nurbs (bnd(ibnd), deriv, deriv2, PTS, 2, rdim);
        end
      end
    end
    
    if (~isempty (boundaries))
      patch_numbers = vertcat (boundaries.patches);
      side_numbers  = vertcat (boundaries.faces);
      if (numel (geometry(1).nurbs.order) > 1)
        for iptc = 1:numel(patch_numbers)
          bnd_nurbs(iptc) = geometry(patch_numbers(iptc)).boundary(side_numbers(iptc)).nurbs;
        end
        boundary_interfaces = nrbmultipatch (bnd_nurbs);
      else
        boundary_interfaces = [];
      end
    end
    
  else
    error('Multiple patches are only implemented for NURBS geometries')
  end
  
end
