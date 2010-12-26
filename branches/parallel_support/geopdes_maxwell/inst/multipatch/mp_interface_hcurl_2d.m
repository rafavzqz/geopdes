% MP_INTERFACE_HCURL_2D: create a global numbering of basis functions in 2d multipatch domains, assigning the correct orientation for tangential traces.
%
%   [glob_num, glob_ndof, dofs_ornt] = mp_interface_hcurl_2d (interfaces, sp);
%
% INPUT:
%
%  interfaces: structure with the information of the interfaces between patches (see mp_geo_read_file)
%  sp:         structure representing the space of discrete functions (see sp_bspline_hcurl_2d_phys)
%
% OUTPUT:
%
%  glob_num:  global numbering of the discrete basis functions
%  glob_ndof: total number of degrees of freedom
%  dofs_ornt: a cell-array with the orientation of the basis functions
%
% Copyright (C) 2010 Rafael Vazquez
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function [glob_num, glob_ndof, dofs_ornt] = mp_interface_hcurl_2d (interfaces, sp)

  dofs_ornt = cell (numel (sp), 1);
  for iptc = 1:numel(sp)
    dofs_ornt{iptc} = ones (1, sp{iptc}.ndof);
  end

  if (~isempty (interfaces))
% Reorder the sides, to avoid conflicts in setting the numbering
    for ind = 1:numel(interfaces)
      if (interfaces(ind).patch1 > interfaces(ind).patch2)
        [interfaces(ind).patch1 interfaces(ind).patch2] = ...
          deal (interfaces(ind).patch2, interfaces(ind).patch1);
        [interfaces(ind).side1 interfaces(ind).side2] = ...
          deal (interfaces(ind).side2, interfaces(ind).side1);
      end
    end

    patches1 = [interfaces.patch1];
    patches2 = [interfaces.patch2];
    sides1   = [interfaces.side1];
    sides2   = [interfaces.side2];
  else
    patches1 = []; sides1 = [];
    patches2 = []; sides2 = [];
  end

  glob_ndof = 0;
  for iptc = 1:numel (sp)
    ind = find (patches2 == iptc);
    new_dofs = setdiff (1:sp{iptc}.ndof, [sp{iptc}.boundary(sides2(ind)).dofs]);
    glob_num{iptc}(new_dofs) = glob_ndof + [1:numel(new_dofs)];
    glob_ndof = glob_ndof + numel(new_dofs);

    for int = find(patches1 == iptc)
      edge = sides1(int);
      intrfc_dofs = sp{iptc}.boundary(edge).dofs;

      nghbr = patches2(int); edge_nghbr = sides2(int);
      nghbr_dofs = sp{nghbr}.boundary(edge_nghbr).dofs;

      if (interfaces(int).ornt == -1)
        nghbr_dofs = fliplr (nghbr_dofs);
        dofs_ornt{nghbr}(nghbr_dofs) = -1;
      end
      glob_num{nghbr}(nghbr_dofs) = glob_num{iptc}(intrfc_dofs);
    end
  end

end
