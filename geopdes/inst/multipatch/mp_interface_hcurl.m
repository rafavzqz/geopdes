% MP_INTERFACE_HCURL: create a global numbering of basis functions in multipatch domains, assigning the correct orientation for tangential traces.
%
%   [glob_num, glob_ndof, dofs_ornt] = mp_interface_hcurl (interfaces, sp);
%
% INPUT:
%
%  interfaces: structure with the information of the interfaces between patches (see mp_geo_read_file)
%  sp:         object representing the space of discrete functions, with the curl-preserving transform (see sp_vector)
%
% OUTPUT:
%
%  glob_num:  global numbering of the discrete basis functions
%  glob_ndof: total number of degrees of freedom
%  dofs_ornt: a cell-array with the orientation of the basis functions
%
% Copyright (C) 2010, 2011, 2015 Rafael Vazquez
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

function [glob_num, glob_ndof, dofs_ornt] = mp_interface_hcurl (interfaces, sp)

  ndim = size (sp{1}.ndof_dir, 2);

  if (~isempty (interfaces))
    glob_num   = cell (numel (sp), 1);
    dofs_ornt  = cell (numel (sp), 1);
% patch_intrfc is only used for ndim == 3
    patch_intrfc = cell (numel (sp), 1);

    ttform = cell (ndim-1, numel (sp), numel (interfaces));
    ppnum  = cell (ndim-1, numel (interfaces));

    for iptc = 1:numel(sp)
      glob_num{iptc} = zeros (1, sp{iptc}.ndof);
      dofs_ornt{iptc} = ones (1, sp{iptc}.ndof);
      patch_intrfc{iptc} = union (find([interfaces.patch1] == iptc), ...
                                  find([interfaces.patch2] == iptc));
    end

    for intrfc = 1:numel(interfaces)
      iptc1  = interfaces(intrfc).patch1;
      iptc2  = interfaces(intrfc).patch2;
      iside1 = interfaces(intrfc).side1;
      iside2 = interfaces(intrfc).side2;
      for idim = 1:ndim-1
        ttform{idim,iptc1, intrfc} = sp{iptc1}.boundary(iside1).comp_dofs{idim};
        nghbr_dofs{idim} = reshape (sp{iptc2}.boundary(iside2).comp_dofs{idim}, ...
				    [sp{iptc2}.boundary(iside2).ndof_dir(idim,:), 1]);
      end
      if (ndim == 2)
        if (interfaces(intrfc).ornt == -1)
          nghbr_dofs{1} = flipud (nghbr_dofs{1});
        end
      elseif (ndim == 3)
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
      end
      for idim = 1:ndim-1
        ttform{idim,iptc2, intrfc} = nghbr_dofs{idim}(:)';
        ppnum{idim,intrfc} = zeros (1, numel(sp{iptc1}.boundary(iside1).comp_dofs{idim}));
      end
    end

    glob_ndof = 0;
% We start with the dofs that do not belong to any interface
    for iptc = 1:numel (sp)
      non_intrfc_dofs = setdiff(1:sp{iptc}.ndof, [ttform{:,iptc,:}]);
      glob_num{iptc}(non_intrfc_dofs) = glob_ndof + (1:numel(non_intrfc_dofs));
      glob_ndof = glob_ndof + numel (non_intrfc_dofs);
    end

% Then we set the interfaces
    for intrfc = 1:numel(interfaces)
      iptc     = interfaces(intrfc).patch1;
      iptc2    = interfaces(intrfc).patch2;

      for idim = 1:ndim-1
        new_dofs{idim} = find (glob_num{iptc}(ttform{idim, iptc, intrfc}) == 0);
        new_dofs_number = glob_ndof + (1:numel (new_dofs{idim}));
        glob_num{iptc}(ttform{idim, iptc, intrfc}(new_dofs{idim})) = new_dofs_number;
        dofs_ornt{iptc}(ttform{idim, iptc, intrfc}(new_dofs{idim})) = 1;
        ppnum{idim,intrfc}(new_dofs{idim}) = new_dofs_number;
        glob_ndof = glob_ndof + numel (new_dofs{idim});
      end

      if (ndim == 2)
        glob_num{iptc2}(ttform{1, iptc2, intrfc}) = ...
             glob_num{iptc}(ttform{1, iptc, intrfc});
        dofs_ornt{iptc2}(ttform{1, iptc2, intrfc}) = ...
             interfaces(intrfc).ornt * dofs_ornt{iptc}(ttform{1, iptc, intrfc});
      elseif (ndim == 3)
        [glob_num, dofs_ornt, ppnum] = set_same_interface (iptc, intrfc, ...
            ttform, new_dofs, interfaces, glob_num, dofs_ornt, ppnum, patch_intrfc);
        [glob_num, dofs_ornt, ppnum] = set_same_patch (iptc, intrfc, ...
            ttform, new_dofs, interfaces, glob_num, dofs_ornt, ppnum, patch_intrfc);

      end
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


function [glob_num, dofs_ornt, ppnum] = set_same_patch (iptc, intrfc, ...
            ttform, new_dofs, interfaces, glob_num, ...
            dofs_ornt, ppnum, patch_intrfc);

  for idim = 1:size(ttform, 1)
    intrfc_dofs{idim} = ttform{idim, iptc, intrfc}(new_dofs{idim});
  end

  for ii = setdiff (patch_intrfc{iptc}, intrfc)
    for idim = 1:size(ttform, 1)
      jdim = mod (idim, 2) + 1;
      [~, pos1{idim}, pos2{idim}] = intersect (ttform{idim, iptc, ii}, intrfc_dofs{idim});
      not_set_p1{idim} = find (ppnum{idim,ii}(pos1{idim}) == 0);

      [~, pos1b{idim}, pos2b{jdim}] = ...
          intersect (ttform{idim, iptc, ii}, intrfc_dofs{jdim});
      not_set_p1b{idim} = find (ppnum{idim, ii}(pos1b{idim}) == 0);
    end

    if (~(all (cellfun (@isempty, not_set_p1))))
      for idim = 1:size(ttform, 1)
        ppnum{idim,ii}(pos1{idim}(not_set_p1{idim})) = ...
          ppnum{idim,intrfc}(new_dofs{idim}(pos2{idim}(not_set_p1{idim})));
        renew_dofs{idim} = pos1{idim}(not_set_p1{idim});
      end
      [glob_num, dofs_ornt, ppnum] = set_same_interface (iptc, ii, ...
          ttform, renew_dofs, interfaces, glob_num, dofs_ornt, ppnum, patch_intrfc);
    elseif (~(all (cellfun (@isempty, not_set_p1b))))
      for idim = 1:size(ttform, 1)
        jdim = mod (idim, 2) + 1;
        ppnum{idim,ii}(pos1b{idim}(not_set_p1b{idim})) = ...
          ppnum{jdim,intrfc}(new_dofs{jdim}(pos2b{jdim}(not_set_p1b{idim})));
        renew_dofs{idim} = pos1b{idim}(not_set_p1b{idim});
      end
      [glob_num, dofs_ornt, ppnum] = set_same_interface (iptc, ii, ...
          ttform, renew_dofs, interfaces, glob_num, dofs_ornt, ppnum, patch_intrfc);
    end
  end
end



function [glob_num, dofs_ornt, ppnum] = set_same_interface (iptc, intrfc, ...
            ttform, new_dofs, interfaces, glob_num, ...
            dofs_ornt, ppnum, patch_intrfc)

  ornt(1) = interfaces(intrfc).ornt1;
  ornt(2) = interfaces(intrfc).ornt2;

  iptc2 = setdiff ([interfaces(intrfc).patch1 interfaces(intrfc).patch2], iptc);
  for idim = 1:size(ttform, 1)
    intrfc_dofs{idim} = ttform{idim, iptc, intrfc}(new_dofs{idim});
    glob_num{iptc2}(ttform{idim, iptc2, intrfc}(new_dofs{idim})) = ...
                              glob_num{iptc}(intrfc_dofs{idim});
    dofs_ornt{iptc2}(ttform{idim, iptc2, intrfc}(new_dofs{idim})) = ...
             ornt(idim) * dofs_ornt{iptc}(intrfc_dofs{idim});
%             interfaces(intrfc).ornt(idim) * dofs_ornt{iptc}(intrfc_dofs{idim});
  end

  [glob_num, dofs_ornt, ppnum] = set_same_patch (iptc2, intrfc, ...
           ttform, new_dofs, interfaces, glob_num, ...
           dofs_ornt, ppnum, patch_intrfc);

end
