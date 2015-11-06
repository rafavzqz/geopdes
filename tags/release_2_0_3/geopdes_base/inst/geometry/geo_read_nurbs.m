% GEO_READ_NURBS: create a geometry structure from a file with a format similar to the one of the NURBS toolbox.
%
% The structure of the geometry file will be available in the documentation
%
% geometry = geo_read_nurbs (filename)
%
% INPUT :
%
%  filename: the name of the file to be read
%
% OUTPUT:
%
%  geometry: a structure containing the following information
%            map:     a function handle to evaluate the parametrization
%            map_der: a function handle to evaluate the derivatives of the parametrization
%            nurbs:   a structure compatible with the NURBS toolbox
%
% Copyright (C) 2009 Carlo de Falco
% Copyright (C) 2010, 2011 Rafael Vazquez
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

function [geom, npatches] = geo_read_nurbs (filename)

fid = fopen (filename, 'r');
if (fid < 0)
  error ('geo_read_nurbs: cannot open file %s', filename);
end

line = fgetl (fid);
while (line ~= -1)
  if (line(1) ~= '#')
    vec = str2num(line);
    dim = vec(1);
    if (numel (vec) > 1)
      npatches = vec(2);
    else
      npatches = 1;
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
  
  if (dim == 2)
    geom(iptc).map      = @(PTS) geo_2d_nurbs (geom(iptc).nurbs, PTS, 0);
    geom(iptc).map_der  = @(PTS) geo_2d_nurbs (geom(iptc).nurbs, PTS, 1);
    geom(iptc).map_der2 = @(PTS) geo_2d_nurbs (geom(iptc).nurbs, PTS, 2);
  elseif (dim == 3)
    geom(iptc).map     = @(PTS) geo_3d_nurbs (geom(iptc).nurbs, PTS, 0);
    geom(iptc).map_der = @(PTS) geo_3d_nurbs (geom(iptc).nurbs, PTS, 1);
  end
end

fclose (fid);

end
