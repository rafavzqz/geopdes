% MP_INTERFACE_HCURL_2D: create a global numbering of basis functions in 2d multipatch domains, assigning the correct orientation for tangential traces.
%
%   [glob_num, glob_ndof, dofs_ornt] = mp_interface_hcurl_2d (interfaces, sp);
%
% INPUT:
%
%  interfaces: structure with the information of the interfaces between patches (see mp_geo_read_file)
%  sp:         object representing the space of discrete functions (see sp_vector_2d_curl_transform)
%
% OUTPUT:
%
%  glob_num:  global numbering of the discrete basis functions
%  glob_ndof: total number of degrees of freedom
%  dofs_ornt: a cell-array with the orientation of the basis functions
%
% Copyright (C) 2010, 2011 Rafael Vazquez
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

  if (~isempty (interfaces))
    glob_ndof  = 0;
    glob_num   = cell (numel (sp), 1);
    ttform     = cell (numel (sp), numel (interfaces));
    dofs_ornt  = cell (numel (sp), 1);
    for iptc = 1:numel(sp)
      glob_num{iptc} = zeros (1, sp{iptc}.ndof);
      dofs_ornt{iptc} = ones (1, sp{iptc}.ndof);
    end

    for intrfc = 1:numel(interfaces)
      iptc1  = interfaces(intrfc).patch1;
      iptc2  = interfaces(intrfc).patch2;
      iside1 = interfaces(intrfc).side1;
      iside2 = interfaces(intrfc).side2;
      ttform{iptc1, intrfc} = sp{iptc1}.boundary(iside1).dofs;
      if (interfaces(intrfc).ornt == 1)
        ttform{iptc2, intrfc} = sp{iptc2}.boundary(iside2).dofs;
      else
        ttform{iptc2, intrfc} = fliplr (sp{iptc2}.boundary(iside2).dofs);
      end
      ppnum{intrfc} = zeros (1, sp{iptc1}.boundary(iside1).ndof);
    end

    glob_ndof = 0;
% We start with the dofs that do not belong to any interface
    for iptc = 1:numel (sp)
      non_intrfc_dofs = setdiff(1:sp{iptc}.ndof, [ttform{iptc,:}]);
      glob_num{iptc}(non_intrfc_dofs) = glob_ndof + (1:numel(non_intrfc_dofs));
      glob_ndof = glob_ndof + numel (non_intrfc_dofs);
    end

% Then we set the interfaces
    for intrfc = 1:numel(interfaces)
      iptc     = interfaces(intrfc).patch1;
      iptc2    = interfaces(intrfc).patch2;
      new_dofs_number = glob_ndof + (1:numel (ttform{iptc, intrfc}));
      glob_num{iptc}(ttform{iptc, intrfc}) = new_dofs_number;
      dofs_ornt{iptc}(ttform{iptc, intrfc}) = 1;

      glob_num{iptc2}(ttform{iptc2, intrfc}) = glob_num{iptc}(ttform{iptc, intrfc});
      dofs_ornt{iptc2}(ttform{iptc2, intrfc}) = ...
               interfaces(intrfc).ornt * dofs_ornt{iptc}(ttform{iptc, intrfc});

      glob_ndof = glob_ndof + numel (new_dofs_number);

    end

  else
    glob_ndof = 0;
    glob_num = cell (numel (sp), 1);
    dofs_ornt  = cell (numel (sp), 1);
    for iptc = 1:numel(sp)
      glob_num{iptc} = glob_ndof + (1:sp{iptc}.ndof);
      glob_ndof = glob_ndof + sp{iptc}.ndof;
      dofs_ornt{iptc} = ones (1, sp{iptc}.ndof);
    end
  end    

end
