% VERTICES_STRUCT: Compute interfaces and vertices struct, necessary to build the C1 basis
%
%     [edges, vertices] = vertices_struct (geometry, interfaces, boundaries, boundary_interfaces)
%
% INPUT:
%
%    geometry:   geometry struct (see mp_geo_load)
%    interfaces: information of connectivity between patches (see mp_geo_load)
%    boundaries: information about the boundary of the domain (see mp_geo_load)
%    boundary_interfaces: information of connectivity between boundary patches (see mp_geo_load)
%
% OUTPUT:
%
%    interfaces: struct array similar to interfaces, but it also contains boundary edges
%    vertices:   struct array. For each vertex it contains:
%      - valence_p: number of patches touching the vertex
%      - valence_e: number of edges touching the vertex
%      - edges:     global numbering of the edges
%      - patches:   global numbering of the patches
%      - patch_reorientation: operations necessary to orient the patch as
%             in the standard configuration (reverse x, reverse y, transpose)
%      - edge_orientation: orientation of the edge in interfaces with
%             respect to the standard configuration in the vertex
%      - boundary_vertex: true if the vertex is a boundary vertex
%
% Copyright (C) 2019-2022 Rafael Vazquez
% Copyright (C) 2019-2021 Cesare Bracco
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

function [interfaces, vertices] = vertices_struct (geometry, interfaces_i, boundaries, boundary_interfaces)

for iptc = 1:numel(geometry)
  msh{iptc} = msh_cartesian ({[0 1] [0 1]}, {0.5 0.5}, {1 1}, geometry(iptc), 'der2', false);
  space{iptc} = sp_bspline ({[0 0 1 1] [0 0 1 1]}, [1 1], msh{iptc});
  scalar_spaces{1} = sp_bspline ({[0 1] [0 0 1 1]}, [0 1], msh{iptc});
  scalar_spaces{2} = sp_bspline ({[0 0 1 1] [0 1]}, [1 0], msh{iptc});
  sp_curl{iptc} = sp_vector (scalar_spaces, msh{iptc}, 'curl-preserving');
end
msh = msh_multipatch (msh, boundaries);
space = sp_multipatch (space, msh, interfaces_i, boundary_interfaces);
sp_curl = sp_multipatch (sp_curl, msh, interfaces_i, boundary_interfaces);

N_int = numel (interfaces_i);
interfaces = interfaces_i;
interfaces_bnd = interfaces_i;
icont = 0;
for hh = 1:numel(boundaries)
  for jj = 1:numel(boundaries(hh).patches)
    icont = icont + 1;
    interfaces(N_int+icont).patch1 = boundaries(hh).patches(jj);
    interfaces(N_int+icont).side1 = boundaries(hh).faces(jj);
    interfaces(N_int+icont).patch2 = [];
    interfaces(N_int+icont).side2 = [];
    interfaces(N_int+icont).ornt = [];    
    interfaces_bnd(N_int+icont).patch1 = boundaries(hh).patches(jj);
    interfaces_bnd(N_int+icont).side1 = boundaries(hh).faces(jj);
    interfaces_bnd(N_int+icont).patch2 = 0;
    interfaces_bnd(N_int+icont).side2 = 0;
  end
end

% Find the operations for reparametrization
if (msh.ndim == 2 && msh.rdim == 2)
  for ii = 1:numel(interfaces)
    interfaces(ii).operations = reorientation_edge_planar (interfaces(ii), geometry);
  end
elseif (msh.ndim == 2 && msh.rdim == 3)
  % In the case of 3D surfaces, check that all patches have the same orientation as the surface
  check_orientation (geometry, interfaces)
  
  for ii = 1:numel(interfaces)
    interfaces(ii).operations = reorientation_edge_3dsurface (interfaces(ii), geometry);
  end
end

% Correspondence between interfaces and edges of the space
int2sp = zeros (sp_curl.ndof, 1);
side2dof = [3 4 1 2]; dof2side = side2dof;
for ii = 1:numel(interfaces)
  patch1 = interfaces(ii).patch1;
  side1 = interfaces(ii).side1;
  int2sp(ii) = sp_curl.gnum{patch1}(side2dof(side1));
end
[~,sp2int] = ismember(1:sp_curl.ndof, int2sp);

% Incidence matrices
C = sparse (sp_curl.ndof, space.ndof);
C_local = [1 0 1 0; -1 0 0 1; 0 1 -1 0; 0 -1 0 -1].';
for iptc = 1:space.npatch
  C(sp_curl.gnum{iptc},space.gnum{iptc}) = spdiags(sp_curl.dofs_ornt{iptc}(:),0,4,4) * C_local;
end

vertices = struct('valence_p', [], 'valence_e', [], 'edges', [], 'patches', [], ...
  'patch_reorientation', [], 'edge_orientation', [], 'boundary_vertex', []);
% Inner and boundary vertex are done in a different way
% For inner vertices, I set the patch with lowest number as the first one
if (isempty (space.boundary))
  boundary_vertices = [];
else
  boundary_vertices = space.boundary.dofs(:).';
end
for ivert = 1:space.ndof
  boundary_vertex = ismember (ivert, boundary_vertices);
  sp_edges = find(C(:,ivert));
  edges = sp2int(sp_edges);
  patches_aux = [interfaces_bnd(edges).patch1; interfaces_bnd(edges).patch2];
  sides_aux = [interfaces_bnd(edges).side1; interfaces_bnd(edges).side2];
  patches = setdiff (unique (patches_aux(:)).', 0);
  position = cellfun (@(gnum) find(gnum==ivert), space.gnum(patches));
  valence = numel (patches);
  gnum_curl = sp2int(cell2mat (sp_curl.gnum(patches)));
  gnum_curl = gnum_curl(:,dof2side);
  
  patch_reorientation = zeros (valence, 3);
  for iptc = 1:valence
    patch_reorientation(iptc,1:3) = reorientation_vertex (position(iptc), geometry(patches(iptc)));
  end

% Determine the first edge and patch to reorder the patches and edges
  if (boundary_vertex)
    [current_edge, current_patch, last_edge] = ...
      find_first_edge_and_patch_boundary (edges, patch_reorientation, sp_curl.boundary.dofs, gnum_curl, sp2int);
    reordering_edges = zeros (1, valence+1);
    reordering_edges(valence+1) = last_edge;
    reordering_patches = zeros (1, valence);
  else
    [current_edge, current_patch] = ...
      find_first_edge_and_patch_interior (patches, edges, patches_aux, sides_aux, patch_reorientation);
    reordering_edges = zeros (1, valence);
    reordering_patches = zeros (1, valence);
  end

% Starting from first edge and patch, reorder the other patches and edges
% reordering_edges and reordering_patches are local indices
  reordering_patches(1) = current_patch;
  reordering_edges(1) = current_edge;
  nn = 1;
  while (nn < valence)
    [~,edges_on_patch] = find (patches_aux == patches(current_patch));
    next_edge = setdiff (edges_on_patch, current_edge);
    [patches_on_edge,~] = find (gnum_curl == edges(next_edge));
    next_patch = setdiff (patches_on_edge, current_patch);
    
    nn = nn + 1;
    reordering_edges(nn) = next_edge;
    reordering_patches(nn) = next_patch;
    current_patch = next_patch;
    current_edge = next_edge;
  end
  edges = edges(reordering_edges);
  patches = patches(reordering_patches);

% Check the orientation of each edge compared to the one in interfaces
  if (boundary_vertex)
    edge_orientation = 2*([interfaces_bnd(edges).patch2] == [patches 0]) - 1;
  else
    edge_orientation = 2*([interfaces(edges).patch2] == patches) - 1;
  end
  
  vertices(ivert).valence_p = valence;
  vertices(ivert).valence_e = numel(edges);
  vertices(ivert).edges = edges;
  vertices(ivert).patches = patches;
  vertices(ivert).patch_reorientation = patch_reorientation(reordering_patches,:);
  vertices(ivert).edge_orientation = edge_orientation;
  vertices(ivert).boundary_vertex = boundary_vertex;
end

end


function check_orientation (geometry, interfaces)
  Nint = find ([interfaces.patch2]); 
  if (~isempty (Nint))
    Nint = Nint(end);
  else
    Nint = 0;
  end

  coords_on_side = {0 0.5; 1 0.5; 0.5 0; 0.5 1}; %4 rows for the sides, 2 columns for the coordinates
  
  for iedge = 1:Nint
    patches = [interfaces(iedge).patch1 interfaces(iedge).patch2];
    nrb_patches = [geometry(patches).nurbs];
    sides = [interfaces(iedge).side1 interfaces(iedge).side2];
    jac_side1 = geometry(patches(1)).map_der(coords_on_side(sides(1),:));
    jac_side2 = geometry(patches(2)).map_der(coords_on_side(sides(2),:));
    normal1 = cross (jac_side1(:,1), jac_side1(:,2)); normal1 = normal1 / norm(normal1);
    normal2 = cross (jac_side2(:,1), jac_side2(:,2)); normal2 = normal2 / norm(normal2);
    if (max (abs (normal1 - normal2)) > 1e-13)
      mssg = sprintf('The patches at interface %d do not follow a global orientation. The involved patches are patch %d (side %d), and patch %d (side %d)', iedge, patches(1), sides(1), patches(2), sides(2));
      error (mssg);
    end
  end
end


function operations = reorientation_edge_planar (interface, geometry)
  patches = [interface.patch1 interface.patch2];
  nrb_patches = [geometry(patches).nurbs];
  sides = [interface.side1 interface.side2];
% Change orientation of first patch
  jac = geometry(patches(1)).map_der({rand(1),rand(1)});
  jacdet = geopdes_det__ (jac);
  if (sides(1) == 2)
    if (jacdet < 0)
      operations(1,:) = [1 0 0];
    else
      operations(1,:) = [1 1 0];
    end
  elseif (sides(1) == 3)
    if (jacdet < 0)
      operations(1,:) = [0 0 1];
    else
      operations(1,:) = [1 0 1];
    end
  elseif (sides(1) == 4)
    if (jacdet < 0)
      operations(1,:) = [1 1 1];
    else
      operations(1,:) = [0 1 1];
    end
  elseif (sides(1) == 1)
    if (jacdet < 0)
      operations(1,:) = [0 1 0];
    else
      operations(1,:) = [0 0 0];
    end
  end

% Change orientation of second patch, only for inner edges
  if (numel(nrb_patches) == 2)
    jac = geometry(patches(2)).map_der({rand(1),rand(1)});
    jacdet = geopdes_det__ (jac);
    if (sides(2) == 1)
      if (jacdet < 0)
        operations(2,:) = [0 0 1];
      else
        operations(2,:) = [0 1 1];
      end
    elseif (sides(2) == 2)
      if (jacdet < 0)
        operations(2,:) = [1 1 1];
      else
        operations(2,:) = [1 0 1];
      end
    elseif (sides(2) == 3)
      if (jacdet < 0)
        operations(2,:) = [1 0 0];
      else
        operations(2,:) = [0 0 0];
      end
    elseif (sides(2) == 4)
      if (jacdet < 0)
        operations(2,:) = [0 1 0];
      else
        operations(2,:) = [1 1 0];
      end
    end
  end
end

function operations = reorientation_edge_3dsurface (interface, geometry)
  patches = [interface.patch1 interface.patch2];
  nrb_patches = [geometry(patches).nurbs];
  sides = [interface.side1 interface.side2];
% Change orientation of first patch
  if (sides(1) == 2)
    operations(1,:) = [1 1 0];
  elseif (sides(1) == 3)
    operations(1,:) = [1 0 1];
  elseif (sides(1) == 4)
    operations(1,:) = [0 1 1];
  elseif (sides(1) == 1)
    operations(1,:) = [0 0 0];
  end

% Change orientation of second patch, only for inner edges
  if (numel(nrb_patches) == 2)
    if (sides(2) == 1)
      operations(2,:) = [0 1 1];
    elseif (sides(2) == 2)
      operations(2,:) = [1 0 1];
    elseif (sides(2) == 3)
      operations(2,:) = [0 0 0];
    elseif (sides(2) == 4)
      operations(2,:) = [1 1 0];
    end
  end
end


function operations = reorientation_vertex (position, geom_patch)

  operations = zeros (1,3);
  switch position
    case 2
      operations(1) = 1;
    case 3
      operations(2) = 1;
    case 4
      operations(1:2) = [1 1];
  end
  
  if (geom_patch.rdim == 2)
    jac = geom_patch.map_der({rand(1),rand(1)});
    jacdet = geopdes_det__ (jac);
    if ((jacdet > 0 && sum(operations) == 1) || (jacdet < 0 && sum(operations) ~= 1))
      operations(3) = 1;
    end
  elseif (geom_patch.rdim == 3)
    if (sum(operations) == 1)
      operations(3) = 1;
    end
  end

end

% Determine the first edge and patch to reorder the patches and edges
% These two functions use local numbering
function [first_edge, first_patch] = find_first_edge_and_patch_interior (patches, edges, patches_aux, sides_aux, patch_reorientation)
  first_patch = 1;
  reornt = patch_reorientation(1,:);
  if (all (reornt == [0 0 0]) || all (reornt == [1 0 0]))
    side = 3;
  elseif (all (reornt == [0 1 0]) || all (reornt == [1 1 0]))
    side = 4;
  elseif (all (reornt == [0 0 1]) || all (reornt == [0 1 1]))
    side = 1;
  else
    side = 2;
  end

  indices = find (patches_aux == patches(first_patch));
  edge_aux = find (sides_aux(indices) == side);
  [~,first_edge] = ind2sub([2,numel(edges)], indices(edge_aux));
end

function [first_edge, first_patch, last_edge] = ...
           find_first_edge_and_patch_boundary (edges, patch_reorientation, curl_bnd_dofs, gnum_curl, sp2int)
  [~,~,bnd_edges] = intersect (sp2int(curl_bnd_dofs), edges);
  if (numel (bnd_edges) ~=2)
    error ('The number of boundary edges near a vertex should be equal to two')
  end
  for ii = 1:2
    [bnd_patches(ii), bnd_sides(ii)] = find (gnum_curl == edges(bnd_edges(ii)));
  end

  reornt = patch_reorientation(bnd_patches(1),:);
  if ((bnd_sides(1) == 3 && (all (reornt == [0 0 0]) || all (reornt == [1 0 0]))) ...
      || (bnd_sides(1) == 4 && (all (reornt == [0 1 0]) || all (reornt == [1 1 0]))) ...
      || (bnd_sides(1) == 1 && (all (reornt == [0 0 1]) || all (reornt == [0 1 1]))) ...
      || (bnd_sides(1) == 2 && (all (reornt == [1 0 1]) || all (reornt == [1 1 1]))))
    first_edge = bnd_edges(1);
    first_patch = bnd_patches(1);
    last_edge = bnd_edges(2);
  else
    first_edge = bnd_edges(2);
    first_patch = bnd_patches(end);
    last_edge = bnd_edges(1);
  end
end
