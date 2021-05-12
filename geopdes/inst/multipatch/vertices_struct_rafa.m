%INPUT: boundaries and interfaces as given by mp_geo_load
%OUTPUT: structures containing 1) for each interfaces (included boundaries) indices of the djacent patches 
%and corresponding incides of sides in the patches; 2)for each vertex the list of interfaces
%containing it and the corresponding index of the point in the interface (1=left/bottom, 2=right/top)

function [interfaces, vertices] = vertices_struct(boundaries, interfaces_i, geometry, boundary_interfaces)

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
for hh = 1:numel(boundaries)
  interfaces(N_int+hh).patch1 = boundaries(hh).patches;
  interfaces(N_int+hh).side1 = boundaries(hh).faces;
  interfaces_bnd(N_int+hh).patch1 = boundaries(hh).patches;
  interfaces_bnd(N_int+hh).side1 = boundaries(hh).faces;
  interfaces_bnd(N_int+hh).patch2 = 0;
  interfaces_bnd(N_int+hh).side2 = 0;
end

% Find the operations for reorientation
for ii = 1:numel(interfaces)
  interfaces(ii).operations = reorientation_edge (interfaces(ii), geometry);
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
boundary_vertices = space.boundary.dofs(:).';
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
    patch_reorientation(iptc,1:3) = reorientation_vertex (position(iptc), geometry(patches(iptc)).nurbs);
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

function operations = reorientation_edge (interface, geometry)
  patches = [interface.patch1 interface.patch2];
  nrb_patches = [geometry(patches).nurbs];
  sides = [interface.side1 interface.side2];
% Change orientation of first patch
  [~,jac] = nrbdeval (nrb_patches(1), nrbderiv(nrb_patches(1)), {rand(1) rand(1)});
  jacdet = jac{1}(1) * jac{2}(2) - jac{1}(2) * jac{2}(1);
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
    [~,jac] = nrbdeval (nrb_patches(2), nrbderiv(nrb_patches(2)), {rand(1) rand(1)});
    jacdet = jac{1}(1) * jac{2}(2) - jac{1}(2) * jac{2}(1);
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


function operations = reorientation_vertex (position, nrb)

  operations = zeros (1,3);
  switch position
    case 2
      operations(1) = 1;
      nrb = nrbreverse (nrb, 1);
    case 3
      operations(2) = 1;
      nrb = nrbreverse (nrb, 2);
    case 4
      operations(1:2) = [1 1];
      nrb = nrbreverse (nrb, [1 2]);
  end
  [~,jac] = nrbdeval (nrb, nrbderiv(nrb), {rand(1) rand(1)});
  jacdet = jac{1}(1) * jac{2}(2) - jac{1}(2) * jac{2}(1);
  
  if (jacdet < 0)
    operations(3) = 1;
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