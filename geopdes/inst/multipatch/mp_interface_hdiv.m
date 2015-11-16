% MP_INTERFACE_HDIV: create a global numbering of basis functions in multipatch domains, assigning the correct orientation for normal traces.
%
%   [glob_num, glob_ndof, dofs_ornt] = mp_interface_hdiv (interfaces, sp, msh);
%
% INPUT:
%
%  interfaces: structure with the information of the interfaces between patches (see mp_geo_read_file)
%  sp:         object representing the space of discrete functions, with a div-preserving transform (see sp_vector)
%
% OUTPUT:
%
%  glob_num:  global numbering of the discrete basis functions
%  glob_ndof: total number of degrees of freedom
%  dofs_ornt: a cell-array with the orientation of the basis functions
%
% Copyright (C) 2010, 2011, 2014, 2015 Rafael Vazquez
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

function [glob_num, glob_ndof, dofs_ornt] = mp_interface_hdiv (interfaces, sp, msh)

  ndim = msh.ndim;

  if (~isempty (interfaces))
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
      ind1 = ceil(iside1/2); %[1 1 2 2 3 3];
      ind2 = ceil(iside2/2); %[1 1 2 2 3 3];
      ttform{iptc1, intrfc} = sp{iptc1}.boundary(iside1).dofs;

      if (ndim == 2)
        if (interfaces(intrfc).ornt == 1)
          ttform{iptc2, intrfc} = sp{iptc2}.boundary(iside2).dofs;
        else
          ttform{iptc2, intrfc} = fliplr (sp{iptc2}.boundary(iside2).dofs);
        end
      elseif (ndim == 3)
        nghbr_dofs = reshape (sp{iptc2}.boundary(iside2).dofs, ...
                      sp{iptc2}.boundary(iside2).ndof_dir);
        ttform{iptc2, intrfc} = sp{iptc2}.boundary(iside2).dofs;
        if (interfaces(intrfc).flag == -1)
          nghbr_dofs = nghbr_dofs';
        end
        if (interfaces(intrfc).ornt1 == -1)
          nghbr_dofs = flipud (nghbr_dofs);
        end
        if (interfaces(intrfc).ornt2 == -1)
          nghbr_dofs = fliplr (nghbr_dofs);
        end
        ttform{iptc2, intrfc} = nghbr_dofs(:)';        
      end
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
      iptc   = interfaces(intrfc).patch1;
      iptc2  = interfaces(intrfc).patch2;
      iside1 = interfaces(intrfc).side1;
      iside2 = interfaces(intrfc).side2;
      new_dofs_number = glob_ndof + (1:numel (ttform{iptc, intrfc}));
      glob_num{iptc}(ttform{iptc, intrfc}) = new_dofs_number;
      dofs_ornt{iptc}(ttform{iptc, intrfc}) = 1;

% Compute the Jacobian for both patches (overkilling)
      msh_el1 = msh_evaluate_element_list (msh.msh_patch{iptc}, 1);
      msh_el2 = msh_evaluate_element_list (msh.msh_patch{iptc2}, 2);
      jac1_sign = sign (geopdes_det__ (msh_el1.geo_map_jac(:,:,1,1)));
      jac2_sign = sign (geopdes_det__ (msh_el2.geo_map_jac(:,:,1,1)));

      if (mod(iside1,2) == mod(iside2,2))
        ornt = -(jac1_sign * jac2_sign);
      else
        ornt =  jac1_sign * jac2_sign;
      end
      glob_num{iptc2}(ttform{iptc2, intrfc}) = glob_num{iptc}(ttform{iptc, intrfc});
      dofs_ornt{iptc2}(ttform{iptc2, intrfc}) = ...
               ornt * dofs_ornt{iptc}(ttform{iptc, intrfc});

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
