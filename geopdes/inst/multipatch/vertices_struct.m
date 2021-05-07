function [interfaces, vertices] = vertices_struct(boundaries, interfaces_i)

%INPUT: boundaries and interfaces as given by mp_geo_load
%OUTPUT: structures containing 1) for each interfaces (included boundaries) indices of the djacent patches 
%and corresponding incides of sides in the patches; 2)for each vertex the list of interfaces
%containing it and the corresponding index of the point in the interface (1=left/bottom, 2=right/top)

N_int = numel (interfaces_i);
interfaces = interfaces_i;
for hh = 1:numel(boundaries)
  interfaces(N_int+hh).patch1 = boundaries(hh).patches;
  interfaces(N_int+hh).side1 = boundaries(hh).faces;
end

%initialize vertices
ivertices=[];

% Find all vertices as intersection of two edges
nv = 0;
for ii = 1:numel(interfaces)
  patches_i = [interfaces(ii).patch1 interfaces(ii).patch2];
  sides_i = [interfaces(ii).side1 interfaces(ii).side2];  
  for jj = ii+1:numel(interfaces)
    patches_j = [interfaces(jj).patch1 interfaces(jj).patch2];
    sides_j = [interfaces(jj).side1 interfaces(jj).side2];
    [cpatch, cpatch_indi, cpatch_indj] = intersect (patches_i, patches_j);
    cside_i = sides_i(cpatch_indi);
    cside_j = sides_j(cpatch_indj);
    
    if (numel(cpatch) > 0) % possible vertex found
      flag = false;
      if (cside_i == 1 && cside_j == 3)
        inds = [1 1]; flag = true;
      elseif (cside_i == 1 && cside_j == 4)
        inds = [2 1]; flag = true;
      elseif (cside_i == 2 && cside_j==3)
        inds = [1 2]; flag = true;
      elseif (cside_i == 2 && cside_j==4)
        inds = [2 2]; flag = true;
      elseif (cside_i == 3 && cside_j == 1)
        inds = [1 1]; flag = true;
      elseif (cside_i == 3 && cside_j == 2)
        inds = [2 1]; flag = true;
      elseif (cside_i == 4 && cside_j == 1)
        inds = [1 2]; flag = true;
      elseif (cside_i == 4 && cside_j==2)
        inds = [2 2]; flag = true;
      end
      if (flag)
        nv = nv+1; %increase number of vertices
        ivertices(nv).interfaces = [ii jj];
        ivertices(nv).ind = inds;
      end
    end
  end
end

%vertices=ivertices;

% Remove redundant vertices
vertices(1) = ivertices(1);
nvert = 1;

for ii = 2:numel(ivertices)
  inter_i = ivertices(ii).interfaces;
  ind_i = ivertices(ii).ind;
  info_i = [inter_i' ind_i'];
  flag_new = 1;
  for jj = 1:numel(vertices)
    inter_j = vertices(jj).interfaces;
    ind_j = vertices(jj).ind;
    info_j = [inter_j' ind_j'];
    if (numel(intersect(info_i,info_j,'rows')) > 0)
      [vertices(jj).interfaces,indu,~] = unique([vertices(jj).interfaces ivertices(ii).interfaces]);
      vertices(jj).ind = [vertices(jj).ind ivertices(ii).ind];
      vertices(jj).ind = vertices(jj).ind(indu);
      flag_new = 0;
      break
    end
  end
  if (flag_new)
    nvert = nvert + 1;
    vertices(nvert) = ivertices(ii);
  end
end   

end
