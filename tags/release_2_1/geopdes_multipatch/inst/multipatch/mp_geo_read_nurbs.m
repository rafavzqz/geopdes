% MP_GEO_READ_NURBS: create a multipatch geometry structure and the interfaces structure from a file in a format similar to the one of the NURBS toolbox.
%
% [geometry, boundaries, interfaces, subdomains] = mp_geo_read_nurbs (filename)
%
% INPUT :
%
%  filename: the name of the file to be read
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
% The interface information is based on the paper
% [1] T. Dokken, E. Quak, V. Skytt, Requirements from Isogeometric Analysis for Changes in Product Design Ontologies. Proceedings of the Focus K3D Conference on Semantic 3D Media and Content, Sophia Antipolis (France), 2010.
%
% An explanation of the file format can be found in the file
%  geopdes_multipatch/doc/geo_specs_mp_v10.txt
%
% Copyright (C) 2009 Carlo de Falco
% Copyright (C) 2010, 2011, 2015 Rafael Vazquez
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

function [geom, boundaries, interfaces, subdomains] = mp_geo_read_nurbs (filename)

  fid = fopen (filename, 'r');
  if (fid < 0)
    error ('mp_geo_read_nurbs: cannot open file %s', filename);
  end
  
  line = fgetl (fid);
  while (line ~= -1)
    if (line(1) ~= '#')
      vec = str2num(line);
      ndim = vec(1);
      rdim = vec(2);
      npatches = vec(3);
      if (numel (vec) == 4)
        nintrfc = vec(4);
        nsubd = 0;
        if (nargout == 4)
          warning ('mp_geo_read_nurbs: no subdomain information is given in the file. The patch number will be used')
          for ii = 1:npatches
            subdomains(ii).name = sprintf ('Patch %i', ii);
            subdomains(ii).patches = ii;
          end
        end
      elseif (numel (vec) == 5)
        nintrfc = vec(4);
        nsubd = vec(5);
      else
        nintrfc = 0;
        interfaces = [];
        nsubd = 0;
        subdomains.name = 'SUBDOMAIN 1';
        subdomains.patches = 1;
      end
      break
    end
    line = fgetl (fid);
  end

  geom(1:npatches) = struct ('map', [], 'map_der', [], 'nurbs', []);
  for iptc = 1:npatches
    geom(iptc).nurbs.form = 'B-NURBS';
    geom(iptc).nurbs.dim = 4;

    line = fgetl (fid);
    while (line(1) == '#')
      line = fgetl (fid);
    end
    if (isempty (str2num (line)))
      geom(iptc).name = line;
      line = fgetl (fid);
      while (line(1) == '#')
        line = fgetl (fid);
      end
    end
    vec = str2num (line);
    geom(iptc).nurbs.order = vec+1;

    line = fgetl (fid);
    while (line(1) == '#')
      line = fgetl (fid);
    end
    vec = str2num (line);
    geom(iptc).nurbs.number = vec;

    for idim = 1:ndim
      line = fgetl (fid);
      while (line(1) == '#')
        line = fgetl (fid);
      end
      geom(iptc).nurbs.knots{idim} = str2num (line);
    end

    for idim = 1:rdim
      line = fgetl (fid); 
      while (line(1) == '#')
        line = fgetl (fid);
      end
      cp_coord = reshape(str2num (line), geom(iptc).nurbs.number);
      geom(iptc).nurbs.coefs(idim,:,:,:) = cp_coord;
    end
    line = fgetl (fid);
    while (line(1) == '#')
      line = fgetl (fid);
    end
    geom(iptc).nurbs.coefs(4,:,:,:) = reshape(str2num (line), geom(iptc).nurbs.number);
          
    geom(iptc).map      = @(PTS) geo_nurbs (geom(iptc).nurbs, PTS, 0);
    geom(iptc).map_der  = @(PTS) geo_nurbs (geom(iptc).nurbs, PTS, 1);
    geom(iptc).map_der2 = @(PTS) geo_nurbs (geom(iptc).nurbs, PTS, 2);
  end

  interfaces = [];
  for intrfc = 1:nintrfc
    line = fgetl (fid);
    while (line(1) == '#')
      line = fgetl (fid);
    end
    interfaces(intrfc).ref = line;
    line = fgetl (fid);
    while (line(1) == '#')
      line = fgetl (fid);
    end
    vec = str2num (line);
    interfaces(intrfc).patch1 = vec(1);
    interfaces(intrfc).side1  = vec(2);
    line = fgetl (fid);
    while (line(1) == '#')
      line = fgetl (fid);
    end
    vec = str2num (line);
    interfaces(intrfc).patch2 = vec(1);
    interfaces(intrfc).side2  = vec(2);
    switch (ndim)
     case 2
       line = fgetl (fid);
       while (line(1) == '#')
         line = fgetl (fid);
       end
       interfaces(intrfc).ornt = str2num (line);
     case 3
       line = fgetl (fid);
       while (line(1) == '#')
         line = fgetl (fid);
       end
       vec = str2num (line);
       interfaces(intrfc).flag  = vec(1); 
       interfaces(intrfc).ornt1 = vec(2);
       interfaces(intrfc).ornt2 = vec(3);
    end
  end

  subdomains = [];
  for isub = 1:nsubd
    line = fgetl (fid);
    while (line(1) == '#')
      line = fgetl (fid);
    end
    subdomains(isub).name = line;
    line = fgetl (fid);
    while (line(1) == '#')
      line = fgetl (fid);
    end
    subdomains(isub).patches = str2num (line);    
  end

  boundaries = [];
  line = fgetl (fid); bnd = 0;
  while (line(1) == '#')
    line = fgetl (fid);
  end
  while (line ~= -1)
    bnd = bnd + 1;
    boundaries(bnd).name = line;
    line = fgetl (fid);
    while (line(1) == '#')
      line = fgetl (fid);
    end
    boundaries(bnd).nsides = str2num (line);
    boundaries(bnd).patches = zeros (boundaries(bnd).nsides, 1); 
    boundaries(bnd).faces = zeros (boundaries(bnd).nsides, 1);
    for isides = 1:boundaries(bnd).nsides
      line = fgetl (fid);
      while (line(1) == '#')
        line = fgetl (fid);
      end
      vec = str2num (line);
      boundaries(bnd).patches(isides) = vec(1); 
      boundaries(bnd).faces(isides) = vec(2);
    end
    line = fgetl (fid);
    while (isempty(line) || line(1) == '#')
      line = fgetl (fid);
    end
  end

  if (isempty(boundaries) && npatches == 1)
    for iface = 1:2*ndim
      boundaries(iface).name = ['BOUNDARY ' num2str(iface)];
      boundaries(iface).nsides = 1;
      boundaries(iface).patches = 1;
      boundaries(iface).faces = iface;
    end
  elseif (isempty (boundaries))
    warning ('The boundary information has not been assigned.\n')
  end

  if (isempty(subdomains))
    subdomains.name = 'WHOLE DOMAIN';
    subdomains.patches = 1:npatches;
  end

  fclose (fid);

end
