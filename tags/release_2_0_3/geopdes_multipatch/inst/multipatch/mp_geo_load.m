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
%      2D Case
%               ref:    reference number of each interface
%               patch1: the index of the first coinciding patch, and the index of its interface boundary edge
%               patch2: the index of the second coinciding patch, and the index of its interface boundary edge
%               ornt:   a flag telling if the direction on patch1 matches the direction on patch2
%      3D Case
%               ref:    reference number of each interface
%               patch1: the index of the first coinciding patch, and the index of its interface boundary face
%               patch2: the index of the second coinciding patch, and the index of its interface boundary face
%               flag:   a flag telling if the u- and v- parameter directions on the first face coincide with the u- and v- directions on the second face
%               ornt1:  a flag telling if the u- direction on the first face matches the direction on patch2
%               ornt2:  a flag telling if the v- direction on the first face matches the direction on patch2
%
%   subdomains: an array of structures that contains the following fields
%               name:    the name of the subdomain
%               patches: indices of the patches belonging to the subdomain
%
%
% This format is based on the paper
% [1] T. Dokken, E. Quak, V. Skytt, Requirements from Isogeometric Analysis for Changes in Product Design Ontologies. Proceedings of the Focus K3D Conference on Semantic 3D Media and Content, Sophia Antipolis (France), 2010.
%
% Copyright (C) 2010, 2011 Carlo de Falco, Rafael Vazquez
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

function [geometry, boundaries, interfaces, subdomains] = mp_geo_load (in)
  
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

      if (dim == 2)
        for iptc = 1:npatch
          geometry(iptc).map      =  ...
            @(PTS) geo_2d_nurbs (geometry(iptc).nurbs, PTS, 0);
          geometry(iptc).map_der  =  ...
            @(PTS) geo_2d_nurbs (geometry(iptc).nurbs, PTS, 1);
          geometry(iptc).map_der2 =  ...
            @(PTS) geo_2d_nurbs (geometry(iptc).nurbs, PTS, 2);
        end
      elseif (dim == 3)
        for iptc = 1:npatch
          geometry(iptc).map      =  ...
            @(PTS) geo_3d_nurbs (geometry(iptc).nurbs, PTS, 0);
          geometry(iptc).map_der  =  ...
            @(PTS) geo_3d_nurbs (geometry(iptc).nurbs, PTS, 1);
        end
      end

    elseif (strcmpi (in(end-3:end), '.txt'))
%% load geometry from a txt file
      [geometry, boundaries, interfaces, subdomains] = mp_geo_read_nurbs (in);
    else
      error ('mp_geo_load: unknown file extension');
    end

  else
    error ('mp_geo_load: wrong input type');
  end

end
