% SP_EVAL_PHYS: Compute the value or derivatives of a function from its degrees of freedom, at a given set of points in the physical domain.
%
%   eu = sp_eval_phys (u, space, geometry, pts, [patch_list], [options]);
%
% INPUT:
%     
%     u:           vector of dof weights
%     space:       object defining the discrete space (see sp_multipatch)
%     geometry:    geometry structure (see mp_geo_load)
%     pts:         array (rdim x npts) with coordinates of points
%     npts:        number of points along each parametric direction
%     patch_list:  patch to which each point begins. If empty, the patch
%                   will be found using nrbinverse
%     options:     cell array with the fields to plot
%                   accepted options are 'value' (default), 'gradient' and 'laplacian'
%
% OUTPUT:
%
%     eu:       cell-array with the fields evaluated at the given points.
%     pts_list: list of points that are on the geometry, for which the
%                values are computed.
%  If there is only one output argument, points not on the geometry get a NaN value.
% 
% Copyright (C) 2023 Rafael Vazquez
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

function eu = sp_eval_phys (u, space, geometry, pts, patch_list, options)

  if (numel (u) ~= space.ndof)
    error ('The number of degrees of freedom of the vector and the space do not match')
  end
  if (numel(geometry) ~= space.npatch)
    error ('The number of patches of the space and geometry do not coincide')
  end    
  if (~isfield (geometry, 'nurbs'))
    error ('Only implemented for NURBS geometries')
  end
  if (nargin > 4 && ~isempty(patch_list) && numel(patch_list) ~= size(pts, 2) )
    error ('The number of patches in the list must be equal to the number of points')
  end

  npatch = numel (geometry);

  if (nargin < 6)
    options = {'value'};
  end
  if (~iscell (options))
    options = {options};
  end
  nopts = numel (options);

  ndim = numel (space.sp_patch{1}.knots);
  rdim = geometry.rdim;
  
  nurbs = geometry.nurbs;
  if (ndim == 1)
    for iptc = 1:npatch
      geometry(iptc).nurbs.knots = {nurbs.knots};
    end
  end

  npts = size (pts, 2);
  pts_param = zeros (ndim, npts); flag = false (npts, 1);
  
% This is necessary because the NURBS toolbox works in the 3D space
  if (size(pts, 1) < 3)
    pts(end+1:3,:) = 0;
  end

  if (nargin < 5 || isempty (patch_list))
    patch_list = zeros (1, npts);
    for ipnt = 1:npts
      for iptc = 1:npatch
        [pts_param(:,ipnt), flag(ipnt)] = nrbinverse (geometry(iptc).nurbs, pts(:,ipnt));
        if (flag(ipnt))
          patch_list(ipnt) = iptc;
          break
        end
      end
    end
  end

  vals = cell (1, npatch);
  pts_on_patch = cell (1, npatch);
  for iptc = 1:npatch
    pts_on_patch{iptc} = find (patch_list == iptc);
    if (~isempty (pts_on_patch{iptc}))
      u_ptc = u(space.gnum{iptc});
      vals{iptc} = sp_eval_phys (u_ptc, space.sp_patch{iptc}, geometry(iptc), pts(:,pts_on_patch{iptc}), options);
      if (nopts == 1)
        vals{iptc} = {vals{iptc}};
      end
    end
  end  

  value = false; grad = false; laplacian = false; hessian = false;
  
  for iopt = 1:nopts
    switch (lower (options{iopt}))
      case 'value'
        eu{iopt} = NaN (1, npts);
        eunum{iopt} = {1};
        eusize{iopt} = npts;
        value = true;

      case 'gradient'
        eu{iopt} = NaN (rdim, npts);
        eunum{iopt} = {1:rdim};
        eusize{iopt} = [rdim, npts];
        grad = true;
        
      case 'laplacian'
        eu{iopt} = NaN (1, npts);
        eunum{iopt} = {1};
        eusize{iopt} = npts;
        laplacian = true;

      case 'hessian'
        eu{iopt} = NaN (rdim, rdim, npts);
        eunum{iopt} = {1:rdim, 1:rdim};
        eusize{iopt} = [rdim, rdim, npts];
        hessian = true;
    end
  end

  for iopt = 1:nopts
    for iptc = 1:npatch
      if (~isempty (pts_on_patch{iptc}))
        eu{iopt}(eunum{iopt}{:},pts_on_patch{iptc}) = vals{iptc}{iopt};
      end
    end
  end
  
  if (nargout == 2)
    pts_list = sort (cell2mat (pts_on_patch));
    for iopt = 1:nopts
      eu{iopt} = eu{iopt}(eunum{iopt}{:}, pts_list);
    end
  end
  
  if (any (isnan(eu{1})))
    warning ('Some points are either not on the geometry, or not on the given patch')
  end
  
  if (nopts == 1)
    eu = eu{1};
  end

end
