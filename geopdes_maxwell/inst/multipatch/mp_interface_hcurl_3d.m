% MP_INTERFACE_HCURL_3D: create a global numbering of basis functions in two-dimensional multipatch geometries.
%
%   [glob_num, glob_ndof, dofs_ornt] = mp_interface_hcurl_3d (interfaces, sp);
%
% INPUT:
%
%  interfaces: structure with the information of the interfaces between patches (see mp_geo_read_file)
%  sp:         object representing the space of discrete functions (see sp_vector_3d_curl_transform)
%
% OUTPUT:
%
%  glob_num:  global numbering of the discrete basis functions
%  glob_ndof: total number of degrees of freedom
%  dofs_ornt: a cell-array with the orientation of the basis functions
%
% Copyright (C) 2011 Rafael Vazquez
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

function [glob_num, glob_ndof, dofs_ornt] = mp_interface_hcurl_3d (interfaces, sp)

  if (~isempty (interfaces))
    glob_ndof = 0;
    glob_num = cell (numel (sp), 1);
    patch_intrfc = cell (numel (sp), 1);
    ttformu  = cell (numel (sp), numel (interfaces));
    ttformv  = cell (numel (sp), numel (interfaces));
    ppnumu   = cell (numel (interfaces), 1);
    ppnumv   = cell (numel (interfaces), 1);
    dofs_ornt = cell (numel (sp), 1);
    for iptc = 1:numel(sp)
      glob_num{iptc} = zeros (1, sp{iptc}.ndof);
      dofs_ornt{iptc} = zeros (1, sp{iptc}.ndof);
      patch_intrfc{iptc} = union (find([interfaces.patch1] == iptc), ...
                                  find([interfaces.patch2] == iptc));
    end

    for intrfc = 1:numel(interfaces)
      iptc1  = interfaces(intrfc).patch1;
      iptc2  = interfaces(intrfc).patch2;
      iside1 = interfaces(intrfc).side1;
      iside2 = interfaces(intrfc).side2;
      ttformu{iptc1, intrfc} = sp{iptc1}.boundary(iside1).comp_dofs{1};
      ttformv{iptc1, intrfc} = sp{iptc1}.boundary(iside1).comp_dofs{2};
      nghbr_dofs{1} = reshape (sp{iptc2}.boundary(iside2).comp_dofs{1}, ...
                    sp{iptc2}.boundary(iside2).ndof_dir(1,:));
      nghbr_dofs{2} = reshape (sp{iptc2}.boundary(iside2).comp_dofs{2}, ...
                    sp{iptc2}.boundary(iside2).ndof_dir(2,:));

      if (interfaces(intrfc).flag == -1)
        aux = nghbr_dofs{1};
        nghbr_dofs{1} = nghbr_dofs{2}';
        nghbr_dofs{2} = aux';
      end
      if (interfaces(intrfc).ornt1 == -1)
        nghbr_dofs{1} = flipud (nghbr_dofs{1});
        nghbr_dofs{2} = flipud (nghbr_dofs{2});
      end
      if (interfaces(intrfc).ornt2 == -1)
        nghbr_dofs{1} = fliplr (nghbr_dofs{1});
        nghbr_dofs{2} = fliplr (nghbr_dofs{2});
      end
      ttformu{iptc2, intrfc} = nghbr_dofs{1}(:)';
      ttformv{iptc2, intrfc} = nghbr_dofs{2}(:)';
      ppnumu{intrfc} = zeros (1, numel(sp{iptc1}.boundary(iside1).comp_dofs{1}));
      ppnumv{intrfc} = zeros (1, numel(sp{iptc1}.boundary(iside1).comp_dofs{2}));
    end

    glob_ndof = 0;
% We start with the dofs that do not belong to any interface
    for iptc = 1:numel (sp)
      non_intrfc_dofs = setdiff(1:sp{iptc}.ndof, [ttformu{iptc,:} ttformv{iptc,:}]);
      glob_num{iptc}(non_intrfc_dofs) = glob_ndof + (1:numel(non_intrfc_dofs));
      dofs_ornt{iptc}(non_intrfc_dofs) = 1;
      glob_ndof = glob_ndof + numel (non_intrfc_dofs);
    end

% Then we set the interfaces
    for intrfc = 1:numel(interfaces)
      iptc     = interfaces(intrfc).patch1;
      new_dofs_u = find (glob_num{iptc}(ttformu{iptc, intrfc}) == 0);
      new_dofs_number = glob_ndof + (1:numel (new_dofs_u));
      glob_num{iptc}(ttformu{iptc, intrfc}(new_dofs_u)) = new_dofs_number;
      glob_ndof = glob_ndof + numel (new_dofs_u);
      ppnumu{intrfc}(new_dofs_u) = new_dofs_number;
      dofs_ornt{iptc}(ttformu{iptc, intrfc}(new_dofs_u)) = 1;

      new_dofs_v = find (glob_num{iptc}(ttformv{iptc, intrfc}) == 0);
      new_dofs_number = glob_ndof + (1:numel (new_dofs_v));
      glob_num{iptc}(ttformv{iptc, intrfc}(new_dofs_v)) = new_dofs_number;
      glob_ndof = glob_ndof + numel (new_dofs_v);
      ppnumv{intrfc}(new_dofs_v) = new_dofs_number;
      dofs_ornt{iptc}(ttformv{iptc, intrfc}(new_dofs_v)) = 1;

      [glob_num, dofs_ornt, ppnumu, ppnumv] = set_same_interface (iptc, intrfc, ...
            ttformu, ttformv, new_dofs_u, new_dofs_v, interfaces, glob_num, ...
            dofs_ornt, ppnumu, ppnumv, patch_intrfc);
      [glob_num, dofs_ornt, ppnumu, ppnumv] = set_same_patch (iptc, intrfc, ...
            ttformu, ttformv, new_dofs_u, new_dofs_v, interfaces, glob_num, ...
            dofs_ornt, ppnumu, ppnumv, patch_intrfc);

    end

  else
    glob_ndof = 0;
    glob_num = cell (numel (sp), 1);
    for iptc = 1:numel(sp)
      glob_num{iptc} = glob_ndof + (1:sp{iptc}.ndof);
      glob_ndof = glob_ndof + sp{iptc}.ndof;
      dofs_ornt{iptc} = ones (1, sp{iptc}.ndof);
    end
  end    

end

function [glob_num, dofs_ornt, ppnumu, ppnumv] = set_same_patch (iptc, intrfc, ...
            ttformu, ttformv, new_dofs_u, new_dofs_v, interfaces, glob_num, ...
            dofs_ornt, ppnumu, ppnumv, patch_intrfc);

  intrfc_dofs_u = ttformu{iptc, intrfc}(new_dofs_u);
  intrfc_dofs_v = ttformv{iptc, intrfc}(new_dofs_v);
  for ii = setdiff (patch_intrfc{iptc}, intrfc)
    [dummy, pos1_u, pos2_u] = intersect (ttformu{iptc, ii}, intrfc_dofs_u);
    [dummy, pos1_v, pos2_v] = intersect (ttformv{iptc, ii}, intrfc_dofs_v);
    not_set_p1u = find (ppnumu{ii}(pos1_u) == 0);
    not_set_p1v = find (ppnumv{ii}(pos1_v) == 0);

    [dummy, pos1_ub, pos2_vb] = intersect (ttformu{iptc, ii}, intrfc_dofs_v);
    [dummy, pos1_vb, pos2_ub] = intersect (ttformv{iptc, ii}, intrfc_dofs_u);
    not_set_p1ub = find (ppnumu{ii}(pos1_ub) == 0);
    not_set_p1vb = find (ppnumv{ii}(pos1_vb) == 0);
    if (~isempty (not_set_p1u) || ~isempty (not_set_p1v))
      ppnumu{ii}(pos1_u(not_set_p1u)) = ppnumu{intrfc}(pos2_u(not_set_p1u));
      ppnumv{ii}(pos1_v(not_set_p1v)) = ppnumv{intrfc}(pos2_v(not_set_p1v));
      [glob_num, dofs_ornt, ppnumu, ppnumv] = set_same_interface (iptc, ii, ...
          ttformu, ttformv, pos1_u(not_set_p1u), pos1_v(not_set_p1v), ...
          interfaces, glob_num, dofs_ornt, ppnumu, ppnumv, patch_intrfc);
    elseif (~isempty (not_set_p1ub) || ~isempty (not_set_p1vb))
      ppnumu{ii}(pos1_ub(not_set_p1ub)) = ppnumv{intrfc}(pos2_vb(not_set_p1ub));
      ppnumv{ii}(pos1_vb(not_set_p1vb)) = ppnumu{intrfc}(pos2_ub(not_set_p1vb));
      [glob_num, dofs_ornt, ppnumu, ppnumv] = set_same_interface (iptc, ii, ...
          ttformu, ttformv, pos1_ub(not_set_p1ub), pos1_vb(not_set_p1vb), ...
          interfaces, glob_num, dofs_ornt, ppnumu, ppnumv, patch_intrfc);
    end
  end
end

function [glob_num, dofs_ornt, ppnumu, ppnumv] = set_same_interface (iptc, intrfc, ...
            ttformu, ttformv, new_dofs_u, new_dofs_v, interfaces, glob_num, ...
            dofs_ornt, ppnumu, ppnumv, patch_intrfc)

  intrfc_dofs_u = ttformu{iptc, intrfc}(new_dofs_u);
  intrfc_dofs_v = ttformv{iptc, intrfc}(new_dofs_v);
  iptc2 = setdiff ([interfaces(intrfc).patch1 interfaces(intrfc).patch2], iptc);
  glob_num{iptc2}(ttformu{iptc2, intrfc}(new_dofs_u)) = ...
                            glob_num{iptc}(intrfc_dofs_u);
  glob_num{iptc2}(ttformv{iptc2, intrfc}(new_dofs_v)) = ...
                            glob_num{iptc}(intrfc_dofs_v);
  dofs_ornt{iptc2}(ttformu{iptc2, intrfc}(new_dofs_u)) = ...
             interfaces(intrfc).ornt1 * dofs_ornt{iptc}(intrfc_dofs_u);
  dofs_ornt{iptc2}(ttformv{iptc2, intrfc}(new_dofs_v)) = ...
             interfaces(intrfc).ornt2 * dofs_ornt{iptc}(intrfc_dofs_v);

  [glob_num, dofs_ornt, ppnumu, ppnumv] = set_same_patch (iptc2, intrfc, ...
           ttformu, ttformv, new_dofs_u, new_dofs_v, interfaces, glob_num, ...
           dofs_ornt, ppnumu, ppnumv, patch_intrfc);

end
