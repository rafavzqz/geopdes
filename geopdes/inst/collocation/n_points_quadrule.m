function rule = n_points_quadrule(nnodes)
  
% some stupid quadrule with 0 or 1 or 2 points in [-1 1]

for idir = 1:numel(nnodes)
    switch nnodes(idir)
        case 0
            rule{idir}(1,:) = [];
            rule{idir}(2,:) = NaN;            
        case 1
            rule{idir}(1,:) = 0;
            rule{idir}(2,:) = NaN;                        
        case 2
            rule{idir}(1,:) = [-1/3 1/3];
            rule{idir}(2,:) = [NaN NaN];                                    
        case 3
            rule{idir}(1,:) = [-0.8 0 0.8];
            rule{idir}(2,:) = [NaN NaN NaN];                                    
    end
    
end
