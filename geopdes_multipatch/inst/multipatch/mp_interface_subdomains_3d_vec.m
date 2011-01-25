% MP_INTERFACE_SUBDOMAINS_3D_VEC: a.
%
%   [glob_num, glob_ndof] = mp_interface_subdomains ();
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

function gn_sub1_to_sub2 = ...
  mp_interface_subdomains_3d_vec (bnd_sub1, bnd_sub2, sp1, sp2, gn1, gn2)

  for ind1 = 1:numel (bnd_sub1)

    ind2 = strcmp ({bnd_sub2.name}, bnd_sub1(ind1).name);
    if (~any (ind2) || isempty (bnd_sub1(ind1).flag))
      continue
    end

    iptc1 = bnd_sub1(ind1).patches;
    face1 = bnd_sub1(ind1).faces;

    iptc2 = bnd_sub2(ind2).patches;
    face2 = bnd_sub2(ind2).faces;

    bnd_dofs1_comp1 = reshape (sp1{iptc1}.boundary(face1).comp_dofs{1}, ...
                  sp1{iptc1}.boundary(face1).ndof_dir(1,:));
    bnd_dofs1_comp2 = reshape (sp1{iptc1}.boundary(face1).comp_dofs{2}, ...
                  sp1{iptc1}.boundary(face1).ndof_dir(2,:));
    bnd_dofs1_comp3 = reshape (sp1{iptc1}.boundary(face1).comp_dofs{3}, ...
                  sp1{iptc1}.boundary(face1).ndof_dir(3,:));

    bnd_dofs2_comp1 = reshape (sp2{iptc2}.boundary(face2).comp_dofs{1}, ...
                  sp2{iptc2}.boundary(face2).ndof_dir(1,:));
    bnd_dofs2_comp2 = reshape (sp2{iptc2}.boundary(face2).comp_dofs{2}, ...
                  sp2{iptc2}.boundary(face2).ndof_dir(2,:));
    bnd_dofs2_comp3 = reshape (sp2{iptc2}.boundary(face2).comp_dofs{3}, ...
                  sp2{iptc2}.boundary(face2).ndof_dir(3,:));

    if (bnd_sub1(ind1).flag == -1)
      bnd_dofs1_comp1 = bnd_dofs1_comp1';
      bnd_dofs1_comp2 = bnd_dofs1_comp2';
      bnd_dofs1_comp3 = bnd_dofs1_comp3';
    end
    if (bnd_sub1(ind1).ornt1 == -1)
      bnd_dofs1_comp1 = flipud (bnd_dofs1_comp1);
      bnd_dofs1_comp2 = flipud (bnd_dofs1_comp2);
      bnd_dofs1_comp3 = flipud (bnd_dofs1_comp3);
    end
    if (bnd_sub1(ind1).ornt2 == -1)
      bnd_dofs1_comp1 = fliplr (bnd_dofs1_comp1);
      bnd_dofs1_comp2 = fliplr (bnd_dofs1_comp2);
      bnd_dofs1_comp3 = fliplr (bnd_dofs1_comp3);
    end
    if (bnd_sub2(ind2).flag == -1)
      bnd_dofs2_comp1 = bnd_dofs2_comp1';
      bnd_dofs2_comp2 = bnd_dofs2_comp2';
      bnd_dofs2_comp3 = bnd_dofs2_comp3';
    end
    if (bnd_sub2(ind2).ornt1 == -1)
      bnd_dofs2_comp1 = flipud (bnd_dofs2_comp1);
      bnd_dofs2_comp2 = flipud (bnd_dofs2_comp2);
      bnd_dofs2_comp3 = flipud (bnd_dofs2_comp3);
    end
    if (bnd_sub2(ind2).ornt2 == -1)
      bnd_dofs2_comp1 = fliplr (bnd_dofs2_comp1);
      bnd_dofs2_comp2 = fliplr (bnd_dofs2_comp2);
      bnd_dofs2_comp3 = fliplr (bnd_dofs2_comp3);
    end

    bnd_dofs1 = [bnd_dofs1_comp1(:); bnd_dofs1_comp2(:); bnd_dofs1_comp3(:)];
    bnd_dofs2 = [bnd_dofs2_comp1(:); bnd_dofs2_comp2(:); bnd_dofs2_comp3(:)];

    gn_sub1_to_sub2{ind1,1} = gn1{iptc1}(bnd_dofs1);
    gn_sub1_to_sub2{ind1,2} = gn2{iptc2}(bnd_dofs2);

  end
end
