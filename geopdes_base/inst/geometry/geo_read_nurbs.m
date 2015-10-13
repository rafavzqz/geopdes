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
%            map:      a function handle to evaluate the parametrization
%            map_der:  a function handle to evaluate the derivatives of the parametrization
%            map_der2: a function handle to evaluate the second derivatives of the parametrization
%            nurbs:    a structure compatible with the NURBS toolbox
%
% Copyright (C) 2009 Carlo de Falco
% Copyright (C) 2010, 2011, 2015 Rafael Vazquez
% Copyright (C) 2014 Elena Bulgarello, Carlo de Falco, Sara Frizziero
% 
% An explanation of the file format can be found in the file
%  geopdes_base/doc/geo_specs_v10.txt
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
    ndim = vec(1);
    if (numel (vec) > 1)
      rdim = vec(2);
    else
      rdim = ndim;
    end
    
    if (numel (vec) > 2)
      npatches = vec(3);
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

  for idim = 1:ndim
    line = fgetl (fid);
    while (line(1) == '#')
      line = fgetl (fid);
    end
    geom(iptc).nurbs.knots{idim} = str2num (line);
  end
  if (ndim == 1)
    geom(iptc).nurbs.knots = geom(iptc).nurbs.knots{1};
  end

  for idim = 1:rdim
    line = fgetl (fid); 
    while (line(1) == '#')
      line = fgetl (fid);
    end
    if (ndim > 1)
      geom(iptc).nurbs.coefs(idim,:,:,:) = reshape (str2num (line), geom(iptc).nurbs.number);
    else
      geom(iptc).nurbs.coefs(idim,:,:,:) = str2num (line);
    end
  end
  line = fgetl (fid);
  while (line(1) == '#')
    line = fgetl (fid);
  end
  if (ndim > 1)
    geom(iptc).nurbs.coefs(4,:,:,:) = reshape (str2num (line), geom(iptc).nurbs.number);
  else
    geom(iptc).nurbs.coefs(4,:,:,:) = str2num (line);
  end
  
  geom(iptc).map      = @(PTS) geo_nurbs (geom(iptc).nurbs, PTS, 0, rdim);
  geom(iptc).map_der  = @(PTS) geo_nurbs (geom(iptc).nurbs, PTS, 1, rdim);
  geom(iptc).map_der2 = @(PTS) geo_nurbs (geom(iptc).nurbs, PTS, 2, rdim);
end

fclose (fid);

end
