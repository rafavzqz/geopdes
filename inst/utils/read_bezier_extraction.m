% READ_BEZIER_EXTRACTION: read a Bezier extraction file, as given by the Autodesk T-splines plug-in for Rhino3d.
%
% USAGE:
%
%   tspline = read_bezier_extraction (filename)
%
% INPUT:
%
%   filename: the name of the file to be read
%
% OUTPUT:
%
%   tspline: structure containing the following fields
%
%     FIELD_NAME     (SIZE)                  DESCRIPTION
%     ndim           (scalar)                dimension of the geometry
%     rdim           (scalar)                dimension of the physical space
%     ndof           (scalar)                number of basis functions / control points
%     nel            (scalar)                number of Bezier elements
%     control_points (4 x ndof vector)       control point coordinates and weights
%     elements       (1 x nel struct)        information for each element (see below)
%     nbnd_sets      (scalar)                number of boundary sets (node and side sets)
%     boundary       (1 x nbnd_sets cell)    boundary information (node and side sets)
%     nelem_sets     (scalar)                number of element sets
%     element_sets   (1 x nelem_sets struct) element_sets information
%     ndom_sets      (scalar)                number of domain sets
%     domain         (1 x ndom_sets struct)  domain information
%
%   The "elements" substructure contains the following fields:
%
%     FIELD_NAME     (SIZE)              DESCRIPTION
%     nsh            (scalar)            number of non-vanishing functions in the element
%     degree         (scalar)            local degree in the element
%     connectivity   (1 x nsh vector)    global numbering of the degrees of freedom
%     extraction     (nsh x nbrn matrix) Bezier extraction operator, with nbrn = prod (degree+1)
%                                          the number of Bernstein polynomials
%
%   For more details, see:
%      M.A. Scott, T.J.R. Hughes, T.W. Sederberg, M.T. Sederberg
%      An integrated approach to engineering design and analysis
%      using the Autodesk T-spline plugin for Rhino3d
%
% Copyright (C) 2012 Rafael Vazquez
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

function tspline = read_bezier_extraction (filename)

fid = fopen (filename, 'r');
if (fid < 0)
  error ('read_bezier_extraction: cannot open file %s', filename);
end

% Initialize.
nbnd_sets = 0; boundary = [];
nelem_sets = 0;
ndom_sets = 0; domain = [];

line = fgetl (fid);
while (all (line ~= -1) || isempty (line))
  if (isempty (line) || line(1) == ' ' || line(1) == '#')
    line = fgetl (fid);
    continue
  else
    spaces = find (isspace (line));
    tag = line(1:spaces(1)-1);
    if (strcmpi (tag, 'type'))
      geo_type = line(spaces(1)+1:end);
      if (strcmpi (geo_type, 'plane'))
        ndim = 2; 
        rdim = 2;
      else
        error ('Only T-spline planar surfaces are implemented in GeoPDEs');
      end
    elseif (strcmpi (tag, 'nodeN'))
      ndof = str2num (line(spaces(1)+1:end));
      nodes = zeros (4, ndof);
    elseif (strcmpi (tag, 'elemN'))
      nel = str2num (line(spaces(1)+1:end));
      nsh_elem = cell (1, nel);
      deg_elem = cell (1, nel);
      ien_elem = cell (1, nel);
      ext_elem = cell (1, nel);
    elseif (strcmpi (tag, 'node'))
      node_read = str2num (line(spaces(1)+1:end));
      nodes(:,1) = node_read';
      for inode = 2:ndof;
        line = fgetl (fid);
        while (line(1) == '#' || line(1) == ' ')
          line = fgetl (fid);
        end
        node_read = str2num (line(6:end));
        nodes(:,inode) = node_read';
      end
      clear node_read
    elseif (strcmpi (tag, 'belem'))
        nsh_elem{1} = str2num (line(spaces(1):spaces(2)));
        deg_elem{1} = str2num (line(spaces(2):end));
        nbernstein = prod (deg_elem{1} + 1);
        line = fgetl (fid);
        while (isempty (line) || line(1) == '#' || line(1) == ' ')
          line = fgetl (fid);
        end
        ien_elem{1} = str2num (line) + 1;
        ext_elem{1} = sparse (nsh_elem{1}, nbernstein);
        for ifun = 1:nsh_elem{1};
          line = fgetl (fid);
          while (isempty (line) || line(1) == '#' || line(1) == ' ')
            line = fgetl (fid);
          end
          ext_elem{1}(ifun,:) = str2num (line);
        end

      for iel = 2:nel
        line = fgetl (fid);
        while (isempty (line) || line(1) == '#' || line(1) == ' ')
          line = fgetl (fid);
        end
        nsh_elem{iel} = str2num (line(spaces(1):spaces(2)));
        deg_elem{iel} = str2num (line(spaces(2):end));
        nbernstein = prod (deg_elem{iel} + 1);
        line = fgetl (fid);
        while (isempty (line) || line(1) == '#' || line(1) == ' ')
          line = fgetl (fid);
        end
        ien_elem{iel} = str2num (line) + 1;
        ext_elem{iel} = sparse (nsh_elem{iel}, nbernstein);
        for ifun = 1:nsh_elem{iel};
          line = fgetl (fid);
          while (isempty (line) || line(1) == '#' || line(1) == ' ')
            line = fgetl (fid);
          end
          ext_elem{iel}(ifun,:) = str2num (line);
        end
      end
      elem = struct ('nsh', nsh_elem, 'degree', deg_elem, ...
		     'connectivity', ien_elem, 'extraction', ext_elem);
      clear nsh_elem deg_elem ien_elem ext_elem
    elseif (strcmpi (tag, 'set'))
      set_type = line(spaces(2)+1:spaces(3)-1);
      if (strcmpi(set_type, 'node'))
        nbnd_sets = nbnd_sets + 1;
        boundary{nbnd_sets}.type = 'node';
        boundary{nbnd_sets}.name = line(spaces(3)+1:spaces(4)-1);
        boundary{nbnd_sets}.dofs = str2num (line(spaces(4):end)) + 1;
      elseif (strcmpi(set_type, 'side'))
        nbnd_sets = nbnd_sets + 1;
        nsides = str2num (line(spaces(1):spaces(2)));
        boundary{nbnd_sets}.type = 'side';
        boundary{nbnd_sets}.name = line(spaces(3)+1:spaces(4)-1);
        boundary{nbnd_sets}.nsides = nsides;
        side_elem = zeros (1, nsides);
        side_pos = cell (1, nsides);
        for iside = 1:nsides
          side_elem(iside) = str2num (line(spaces(2+2*iside):spaces(3+2*iside))) + 1;
% I save only the first character (L(eft), R(ight), B(ottom), T(op))
          side_pos{iside} = line(spaces(3+2*iside)+1);
        end
        boundary{nbnd_sets}.elem = side_elem;
        boundary{nbnd_sets}.position = side_pos;
      elseif (strcmpi(set_type, 'elem'))
        nelem_sets = nelem_sets + 1;
        element_sets(nelem_sets).type = 'elem';
        element_sets(nelem_sets).name = line(spaces(3)+1:spaces(4)-1);
        element_sets(nelem_sets).element_numbers = str2num (line(spaces(4):end)) + 1;
      end
    elseif (strcmpi (tag, 'domain'))
      ndom_sets = ndom_sets + 1;
      domain(ndom_sets).nel = str2num (line(spaces(1):spaces(2)));
      domain(ndom_sets).ndof = str2num (line(spaces(2):end));
      line = fgetl (fid);
      while (isempty (line) || line(1) == '#' || line(1) == ' ')
        line = fgetl (fid);
      end
      domain(ndom_sets).element_list = str2num (line) + 1;
      line = fgetl (fid);
      while (isempty (line) || line(1) == '#' || line(1) == ' ')
        line = fgetl (fid);
      end
      domain(ndom_sets).dofs = str2num (line) + 1;
      for iel = 1:domain(ndom_sets).nel
        line = fgetl (fid);
        while (isempty(line) || line(1) == '#' || line(1) == ' ')
          line = fgetl (fid);
        end
        domain(ndom_sets).connectivity{iel} = str2num (line) + 1;
      end
    end
  end
  line = fgetl (fid);
end

fclose (fid);

% Save all the information in the tspline structure
tspline = struct ('ndim', ndim, 'rdim', rdim, 'ndof', ndof, 'nel', nel, ...
    'control_points', nodes, 'elements', elem);
    
if (nbnd_sets > 0)
  tspline.nbnd_sets = nbnd_sets;
  tspline.boundary = boundary;  
end
if (nelem_sets > 0)
  tspline.nelem_sets = nelem_sets;
  tspline.element_sets = element_sets;
end
if (ndom_sets > 0)
  tspline.ndom_sets = ndom_sets;
  tspline.domain = domain;
end

end
