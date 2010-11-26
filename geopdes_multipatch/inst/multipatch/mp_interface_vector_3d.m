% MP_INTERFACE_VECTOR_3D: create a global numbering of vectorial basis functions in three-dimensional multipatch geometries.
%
%   [glob_num, glob_ndof] = mp_interface_vector_3d (interfaces, sp);
%
% INPUT:
%
%  interfaces: structure with the information of the interfaces between patches (see mp_geo_read_file)
%  sp:         structure representing the space of discrete functions (see sp_bspline_3d_phys)
%
% OUTPUT:
%
%  glob_num:  global numbering of the discrete basis functions
%  glob_ndof: total number of degrees of freedom
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
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

function [glob_num, glob_ndof] = mp_interface_vector_3d (interfaces, sp)

  if (~isempty (interfaces))
% Reorder the sides, to avoid conflicts in setting the numbering
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
  for iptc = 1:numel (sp)
    ind = find (patches2 == iptc);
    new_dofs = setdiff (1:sp{iptc}.ndof, [sp{iptc}.boundary(sides2(ind)).dofs]);
    glob_num{iptc}(new_dofs) = glob_ndof + [1:numel(new_dofs)];
    glob_ndof = glob_ndof + numel(new_dofs);

    for int = find (patches1 == iptc)
      face = sides1(int);
      intrfc_dofs = sp{iptc}.boundary(face).dofs;

      nghbr = patches2(int); face_nghbr = sides2(int);
      nghbr_dofs_comp1 = reshape (sp{nghbr}.boundary(face_nghbr).comp_dofs{1}, ...
                    sp{nghbr}.boundary(face_nghbr).ndof_dir(1,:));
      nghbr_dofs_comp2 = reshape (sp{nghbr}.boundary(face_nghbr).comp_dofs{2}, ...
                    sp{nghbr}.boundary(face_nghbr).ndof_dir(2,:));
      nghbr_dofs_comp3 = reshape (sp{nghbr}.boundary(face_nghbr).comp_dofs{3}, ...
                    sp{nghbr}.boundary(face_nghbr).ndof_dir(3,:));

      if (interfaces(int).flag == -1)
        nghbr_dofs_comp1 = nghbr_dofs_comp1';
        nghbr_dofs_comp2 = nghbr_dofs_comp2';
        nghbr_dofs_comp3 = nghbr_dofs_comp3';
      end
      if (interfaces(int).ornt1 == -1)
        nghbr_dofs_comp1 = flipud (nghbr_dofs_comp1);
        nghbr_dofs_comp2 = flipud (nghbr_dofs_comp2);
        nghbr_dofs_comp3 = flipud (nghbr_dofs_comp3);
      end
      if (interfaces(int).ornt2 == -1)
        nghbr_dofs_comp1 = fliplr (nghbr_dofs_comp1);
        nghbr_dofs_comp2 = fliplr (nghbr_dofs_comp2);
        nghbr_dofs_comp3 = fliplr (nghbr_dofs_comp3);
      end
      nghbr_dofs = [nghbr_dofs_comp1(:) nghbr_dofs_comp2(:) nghbr_dofs_comp3(:)];

      glob_num{nghbr}(nghbr_dofs) = glob_num{iptc}(intrfc_dofs(:));
    end
  end

end
