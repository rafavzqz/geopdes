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
%     patch_list:  patch to which each point begins. If empty, the patch
%                   will be found using nrbinverse
%     options:     cell array with the fields to plot
%                   accepted options for scalars are 'value' (default), 'gradient', 'laplacian', 'hessian'
%                   accepted options for vectors are 'value' (default), 'gradient', 'curl', 'divergence'
%
% OUTPUT:
%
%     eu:       cell-array with the fields evaluated at the given points.
%
%  Points not on the geometry get a NaN value.
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

  is_scalar = isa (space.sp_patch{1}, 'sp_scalar');
  if (is_scalar) 
    ndim = numel (space.sp_patch{1}.knots);
  else
    ndim = numel (space.sp_patch{1}.scalar_spaces{1}.knots);
  end
  rdim = geometry(1).rdim;
  
  if (ndim == 1)
    for iptc = 1:npatch
      geometry(iptc).nurbs.knots = {geometry(iptc).nurbs.knots};
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
      if (~isempty (space.dofs_ornt))
        u_ptc = u_ptc .* space.dofs_ornt{iptc}(:);
      end
      vals{iptc} = sp_eval_phys (u_ptc, space.sp_patch{iptc}, geometry(iptc), pts(:,pts_on_patch{iptc}), options);
      if (nopts == 1)
        vals{iptc} = {vals{iptc}};
      end
    end
  end  

  [eu, eunum] = set_output_sizes (ndim, rdim, npts, space.ncomp, is_scalar, options);

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

function [eu, eunum] = set_output_sizes (ndim, rdim, npts, ncomp, is_scalar, options)

  nopts = numel(options);
  eu = cell (nopts, 1); eunum = eu;

  if (is_scalar)
    for iopt = 1:nopts
      switch (lower (options{iopt}))
        case 'value'
          eu{iopt} = NaN (1, npts);
          eunum{iopt} = {1};
        case 'gradient'
          eu{iopt} = NaN (rdim, npts);
          eunum{iopt} = {1:rdim};
        case 'laplacian'
          eu{iopt} = NaN (1, npts);
          eunum{iopt} = {1};
        case 'hessian'
          eu{iopt} = NaN (rdim, rdim, npts);
          eunum{iopt} = {1:rdim, 1:rdim};
      end
    end
  else
    for iopt = 1:numel(options)
      switch (lower (options{iopt}))
        case 'value'
          eu{iopt} = NaN (ncomp, npts);
          eunum{iopt} = {1:ncomp};
        case 'gradient'
          eu{iopt} = NaN (ncomp, rdim, npts);
          eunum{iopt} = {1:ncomp, 1:rdim};
        case 'curl'
          if (ndim == 2 && rdim == 2)
            eu{iopt} = NaN (1, npts);
            eunum{iopt} = {1};
          elseif (ndim == 3 && rdim == 3)
            eu{iopt} = NaN (ncomp, rdim, npts);
            eunum{iopt} = {1:ncomp};
          end
        case 'divergence'
          eu{iopt} = NaN (1, npts);
          eunum{iopt} = {1};
      end
    end
  end
end
