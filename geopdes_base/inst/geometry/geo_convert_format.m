% GEO_CONVERT_FORMAT: convert a GeoPDEs geometry file from version 0,7 (or 0,6)
%            to version 2,1
%
% geo_convert_format (input_file, output_file)
%
% INPUT :
%
%  input_file:  the name of the file in the old version, to be read
%  output_file: the name of the file in the new version, to be written
%
% Copyright (C) 2009 Carlo de Falco
% Copyright (C) 2010, 2011, 2015 Rafael Vazquez
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

function geo_convert_format (input_file, output_file)

[geom, boundaries, interfaces, subdomains] = mp_geo_read_nurbs_07 (input_file);
nrb = [geom.nurbs];
npatches = numel (nrb);
ndim = numel (nrb(1).order);
if (npatches == 1 && numel (boundaries) == 2*ndim)
  boundaries = [];
end
nrbexport (nrb, interfaces, boundaries, subdomains, output_file)

end


function [geom, boundaries, interfaces, subdomains] = mp_geo_read_nurbs_07 (filename)

  fid = fopen (filename, 'r');
  if (fid < 0)
    error ('geo_convert_format: cannot open file %s', filename);
  end
  
  line = fgetl (fid);
  while (line ~= -1)
    if (line(1) ~= '#')
      vec = str2num(line);
      dim = vec(1);
      npatches = vec(2);
      if (numel (vec) == 3)
        nintrfc = vec(3);
        nsubd = 0;
        if (nargout == 4)
          warning ('geo_convert_format: no subdomain information is given in the file. The patch number will be used')
          for ii = 1:npatches
            subdomains(ii).name = sprintf ('Patch %i', ii);
            subdomains(ii).patches = ii;
          end
        end
      elseif (numel (vec) == 4)
        nintrfc = vec(3);
        nsubd = vec(4);
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

    for idim = 1:dim
      line = fgetl (fid);
      while (line(1) == '#')
        line = fgetl (fid);
      end
      geom(iptc).nurbs.knots{idim} = str2num (line);
    end

    switch (dim)
     case 2
      line = fgetl (fid); 
      while (line(1) == '#')
        line = fgetl (fid);
      end
      cp_x = reshape(str2num (line), geom(iptc).nurbs.number);
      
      line = fgetl (fid);
      while (line(1) == '#')
        line = fgetl (fid);
      end
      cp_y = reshape(str2num (line), geom(iptc).nurbs.number);

      cp_z = zeros(geom(iptc).nurbs.number);
      
      line = fgetl (fid);
      while (line(1) == '#')
        line = fgetl (fid);
      end
      weights = reshape(str2num (line), geom(iptc).nurbs.number);

      geom(iptc).nurbs.coefs(1,:,:) = cp_x;
      geom(iptc).nurbs.coefs(2,:,:) = cp_y;
      geom(iptc).nurbs.coefs(3,:,:) = cp_z;
      geom(iptc).nurbs.coefs(4,:,:) = weights;

     case 3
      line = fgetl (fid);
      while (line(1) == '#')
        line = fgetl (fid);
      end
      cp_x = reshape(str2num (line), geom(iptc).nurbs.number);

      line = fgetl (fid);
      while (line(1) == '#')
        line = fgetl (fid);
      end
      cp_y = reshape(str2num (line), geom(iptc).nurbs.number);

      line = fgetl (fid);
      while (line(1) == '#')
        line = fgetl (fid);
      end
      cp_z = reshape(str2num (line), geom(iptc).nurbs.number);

      line = fgetl (fid);
      while (line(1) == '#')
        line = fgetl (fid);
      end
      weights = reshape(str2num (line), geom(iptc).nurbs.number);

      geom(iptc).nurbs.coefs(1,:,:,:) = cp_x;
      geom(iptc).nurbs.coefs(2,:,:,:) = cp_y;
      geom(iptc).nurbs.coefs(3,:,:,:) = cp_z;
      geom(iptc).nurbs.coefs(4,:,:,:) = weights;

    end
      
    geom(iptc).map      = @(PTS) geo_nurbs (geom(iptc).nurbs, PTS, 0);
    geom(iptc).map_der  = @(PTS) geo_nurbs (geom(iptc).nurbs, PTS, 1);
    geom(iptc).map_der2 = @(PTS) geo_nurbs (geom(iptc).nurbs, PTS, 2);
  end
  
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
    switch (dim)
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
    for iface = 1:2*dim
      boundaries(iface).name = ['BOUNDARY ' num2str(iface)];
      boundaries(iface).nsides = 1;
      boundaries(iface).patches = 1;
      boundaries(iface).faces = iface;
    end
  elseif (isempty (boundaries))
    warning ('The boundary information has not been assigned.\n')
  end

  fclose (fid);

end
