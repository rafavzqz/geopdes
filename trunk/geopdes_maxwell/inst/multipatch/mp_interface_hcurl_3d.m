% MP_INTERFACE_HCURL_3D: create a global numbering of basis functions in 3d multipatch domains, assigning the correct orientation for tangential traces.
%
%   [glob_num, glob_ndof, dofs_ornt] = mp_interface_hcurl_3d (interfaces, sp);
%
% INPUT:
%
%  interfaces: structure with the information of the interfaces between patches (see mp_geo_read_file)
%  sp:         structure representing the space of discrete functions (see sp_bspline_hcurl_3d_phys)
%
% OUTPUT:
%
%  glob_num:   global numbering of the discrete basis functions
%  glob_ndof:  total number of degrees of freedom
%  dofs_ornt:  a cell-array with the orientation of the basis functions
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

function [glob_num, glob_ndof, dofs_ornt] = mp_interface_hcurl_3d (interfaces, sp)

  edges_in_face = [4 12 5 8; 2 10 6 7; 1 9 5 6; 3 11 8 7; 1 3 4 2; 9 11 12 10];

  dofs_ornt = cell (numel (sp), 1);
  for iptc = 1:numel(sp)
    dofs_ornt{iptc} = zeros (1, sp{iptc}.ndof);
  end

  if (~isempty (interfaces))
% Reorder the sides, to avoid conflicts in setting the numbering and orientation
    for ind = 1:numel(interfaces)
      if (interfaces(ind).patch1 > interfaces(ind).patch2)
        [interfaces(ind).patch1 interfaces(ind).patch2] = ...
          deal (interfaces(ind).patch2, interfaces(ind).patch1);
        [interfaces(ind).side1 interfaces(ind).side2] = ...
          deal (interfaces(ind).side2, interfaces(ind).side1);
        if (interfaces(ind).flag == -1)
          [interfaces(ind).ornt1 interfaces(ind).ornt2] = ...
            deal (interfaces(ind).ornt2, interfaces(ind).ornt1);
        end
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
  nedges = 0;
  edge_ornt = zeros (12, numel(sp));

  for iptc = 1:numel(sp)
    ind = find (patches2 == iptc);
    new_dofs = setdiff (1:sp{iptc}.ndof, [sp{iptc}.boundary(sides2(ind)).dofs]);
    glob_num{iptc}(new_dofs) = glob_ndof + [1:numel(new_dofs)];
    glob_ndof = glob_ndof + numel(new_dofs);
    dofs_ornt{iptc}(new_dofs) = 1;

    cntd_edges = edges_in_face(sides2(ind),:);
    new_edges  = setdiff (1:12, cntd_edges(:));
    edge_num{iptc}(new_edges) = nedges + [1:numel(new_edges)];
    nedges = nedges + numel(new_edges);
    edge_ornt(new_edges, iptc) = 1;

    for int = find (patches1 == iptc)
      face = sides1(int);
      intrfc_dofs = sp{iptc}.boundary(face).dofs;

      nghbr = patches2(int); face_nghbr = sides2(int);

      nghbr_dofs{1} = reshape (sp{nghbr}.boundary(face_nghbr).comp_dofs{1}, ...
                    sp{nghbr}.boundary(face_nghbr).ndof_dir(1,:));
      nghbr_dofs{2} = reshape (sp{nghbr}.boundary(face_nghbr).comp_dofs{2}, ...
                    sp{nghbr}.boundary(face_nghbr).ndof_dir(2,:));

      if (interfaces(int).flag == -1)
        aux = nghbr_dofs{1};
        nghbr_dofs{1} = nghbr_dofs{2}';
        nghbr_dofs{2} = aux';
        edg_u = [3 4];
        edg_v = [1 2];
      else
        edg_u = [1 2];
        edg_v = [3 4];
      end

      if (interfaces(int).ornt1 == -1)
        nghbr_dofs{1} = flipud (nghbr_dofs{1});
        nghbr_dofs{2} = flipud (nghbr_dofs{2});
        aux_dofs = nghbr_dofs{1}(:,2:end-1);
        dofs_ornt{nghbr}(aux_dofs) = -1;
      else
        aux_dofs = nghbr_dofs{1}(:,2:end-1);
        dofs_ornt{nghbr}(aux_dofs) = 1;
      end
      if (edge_ornt(edges_in_face(face_nghbr, edg_u(1)), nghbr) == 0)
        edge_ornt(edges_in_face(face_nghbr, edg_u(1)), nghbr) = ...
                                                 interfaces(int).ornt1;
        aux_dofs = nghbr_dofs{1}(:,1);
        dofs_ornt{nghbr}(aux_dofs) = interfaces(int).ornt1;
      end
      if (edge_ornt(edges_in_face(face_nghbr, edg_u(2)), nghbr) == 0)
        edge_ornt(edges_in_face(face_nghbr, edg_u(2)), nghbr) = ...
                                                 interfaces(int).ornt1;
        aux_dofs = nghbr_dofs{1}(:,end);
        dofs_ornt{nghbr}(aux_dofs) = interfaces(int).ornt1;
      end

      if (interfaces(int).ornt2 == -1)
        nghbr_dofs{1} = fliplr (nghbr_dofs{1});
        nghbr_dofs{2} = fliplr (nghbr_dofs{2});
        aux_dofs = nghbr_dofs{2}(2:end-1,:);
        dofs_ornt{nghbr}(aux_dofs) = -1;
      else
        aux_dofs = nghbr_dofs{2}(2:end-1,:);
        dofs_ornt{nghbr}(aux_dofs) = 1;
      end
      if (edge_ornt(edges_in_face(face_nghbr, edg_v(1)), nghbr) == 0)
        edge_ornt(edges_in_face(face_nghbr, edg_v(1)), nghbr) = ...
                                                 interfaces(int).ornt2;
        aux_dofs = nghbr_dofs{2}(1,:);
        dofs_ornt{nghbr}(aux_dofs) = interfaces(int).ornt2;
      end
      if (edge_ornt(edges_in_face(face_nghbr, edg_v(2)), nghbr) == 0)
        edge_ornt(edges_in_face(face_nghbr, edg_v(2)), nghbr) = ...
                                                 interfaces(int).ornt2;
        aux_dofs = nghbr_dofs{2}(end,:);
        dofs_ornt{nghbr}(aux_dofs) = interfaces(int).ornt2;
      end

      nghbr_dofs = [nghbr_dofs{1}(:); nghbr_dofs{2}(:)];

      glob_num{nghbr}(nghbr_dofs) = glob_num{iptc}(intrfc_dofs(:));

      clear nghbr_dofs
    end
  end
end
