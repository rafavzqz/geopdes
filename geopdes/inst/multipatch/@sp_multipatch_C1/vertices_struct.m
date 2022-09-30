function [interfaces, vertices] = vertices_struct(boundaries, interfaces_i)

%INPUT: boundaries and (interior) interfaces as given by mp_geo_load
%OUTPUT: structures containing 1) for each interface (included boundaries) indices of the adjacent patches 
%and corresponding incides of sides in the patches; 2)for each vertex the list of interfaces
%containing it and the corresponding index of the point in the interface (1=left/bottom, 2=right/top)

N_int=numel(interfaces_i);
interfaces=interfaces_i;
for h=1:numel(boundaries)
    interfaces(N_int+h).patch1=boundaries(h).patches;
    interfaces(N_int+h).patch2=[];
    interfaces(N_int+h).side1=boundaries(h).faces;
    interfaces(N_int+h).side2=[];
    interfaces(N_int+h).ornt=[];
end

%initialize vertices
ivertices=[];

nv=0;
for i=1:numel(interfaces)
    patches_i=[interfaces(i).patch1 interfaces(i).patch2];
    sides_i=[interfaces(i).side1 interfaces(i).side2];  
    for j=i+1:numel(interfaces)
        patches_j=[interfaces(j).patch1 interfaces(j).patch2];
        sides_j=[interfaces(j).side1 interfaces(j).side2];
        [cpatch,cpatch_indi,cpatch_indj] = intersect(patches_i,patches_j);
        cside_i=sides_i(cpatch_indi);
        cside_j=sides_j(cpatch_indj);
        if numel(cpatch)>0 %if vertex found
            switch cside_i
                case 1
                    if cside_j==3
                        nv=nv+1; %increase number of vertices
                        ivertices(nv).interfaces=[i j];
                        ivertices(nv).ind=[1 1];
                    elseif cside_j==4
                        nv=nv+1; %increase number of vertices
                        ivertices(nv).interfaces=[i j];
                        ivertices(nv).ind=[2 1];                        
                    end
                case 2
                    if cside_j==3
                        nv=nv+1; %increase number of vertices
                        ivertices(nv).interfaces=[i j];
                        ivertices(nv).ind=[1 2];
                    elseif cside_j==4
                        nv=nv+1; %increase number of vertices
                        ivertices(nv).interfaces=[i j];
                        ivertices(nv).ind=[2 2];                        
                    end                    
                case 3
                    if cside_j==1
                        nv=nv+1; %increase number of vertices
                        ivertices(nv).interfaces=[i j];
                        ivertices(nv).ind=[1 1];
                    elseif cside_j==2
                        nv=nv+1; %increase number of vertices
                        ivertices(nv).interfaces=[i j];
                        ivertices(nv).ind=[2 1];                        
                    end                    
                case 4     
                    if cside_j==1
                        nv=nv+1; %increase number of vertices
                        ivertices(nv).interfaces=[i j];
                        ivertices(nv).ind=[1 2];
                    elseif cside_j==2
                        nv=nv+1; %increase number of vertices
                        ivertices(nv).interfaces=[i j];
                        ivertices(nv).ind=[2 2];                        
                    end                    
            end
        end
    end
end

%vertices=ivertices;

vertices(1)=ivertices(1);
nv=1;

for i=2:numel(ivertices)
inter_i=ivertices(i).interfaces;
ind_i=ivertices(i).ind;     
info_i=[inter_i' ind_i'];
    flag_new=1;
    for j=1:numel(vertices)
        inter_j=vertices(j).interfaces;
        ind_j=vertices(j).ind;
        info_j=[inter_j' ind_j'];
        if numel(intersect(info_i,info_j,'rows'))>0
            [vertices(j).interfaces,indu,~]=unique([vertices(j).interfaces ivertices(i).interfaces]);
            vertices(j).ind=[vertices(j).ind ivertices(i).ind];
            vertices(j).ind=vertices(j).ind(indu);
            flag_new=0;
            break
        end
    end
    if flag_new==1
        nv=nv+1;
        vertices(nv)=ivertices(i);
    end
end   
    

end

