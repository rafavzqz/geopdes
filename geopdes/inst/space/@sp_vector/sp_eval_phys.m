% SP_EVAL_PHYS: Compute the value or derivatives of a function from its degrees of freedom, at a given set of points in the physical domain.
%
%   eu = sp_eval_phys (u, space, geometry, pts, [options]);
%   [eu, pts_list] = sp_eval_phys (u, space, geometry, pts, [options]);
%
% INPUT:
%     
%     u:           vector of dof weights
%     space:       object defining the discrete space (see sp_vector)
%     geometry:    geometry structure (see geo_load)
%     pts:         array (rdim x npts) with coordinates of points
%     options:     cell array with the fields to plot
%                   accepted options are 'value' (default), 'gradient', 'curl' and 'divergence'
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

function [eu, pts_list] = sp_eval_phys (u, space, geometry, pts, options)

  if (numel (u) ~= space.ndof)
    error ('The number of degrees of freedom of the vector and the space do not match')
  end
  if (~isfield (geometry, 'nurbs'))
    error ('Only implemented for NURBS geometries')
  end

  if (nargin < 5)
    options = {'value'};
  end
  if (~iscell (options))
    options = {options};
  end
  nopts = numel (options);

  ndim = numel (space.scalar_spaces{1}.knots);
  rdim = geometry.rdim;
  
  nurbs = geometry.nurbs;

  npts = size (pts, 2);
  pts_param = zeros (ndim, npts); flag = false (npts, 1);
% This is necessary because the NURBS toolbox works in the 3D space
  if (size(pts, 1) < 3)
    pts(end+1:3,:) = 0;
  end
  
  for ipt = 1:npts
    [pts_param(:,ipt), flag(ipt)] = nrbinverse (nurbs, pts(:,ipt));
  end
  
  if (~all(flag))
    warning ('Some of the points are not on the geometry')
    pts_param = pts_param(:,flag);
    npts_on_F = size (pts_param, 2);
    if (nargout == 2)
      npts_out = npts_on_F;
      pts_list = find (flag(:)');
      inds_on_F = 1:npts_out;
    else
      npts_out = npts;
      inds_on_F = find (flag);
    end
  else
    npts_on_F = npts;
    npts_out = npts;
    inds_on_F = 1:npts;
    pts_list = 1:npts;
  end
  
  for idim = 1:ndim
    zeta{idim} = unique (space.scalar_spaces{1}.knots{idim});
    ind{idim} = findspan (numel(zeta{idim})-2, 0, pts_param(idim,:), zeta{idim}) + 1;
  end
      
  value = false; grad = false; curl = false; divergence = false;
  
  for iopt = 1:nopts
    switch (lower (options{iopt}))
      case 'value'
        eu{iopt} = NaN (space.ncomp, npts_out);
        eunum{iopt} = {1:space.ncomp};
        eusize{iopt} = [space.ncomp, npts_out];
        value = true;

      case 'gradient'
        eu{iopt} = NaN (space.ncomp, rdim, npts_out);
        eunum{iopt} = {1:space.ncomp, 1:rdim};
        eusize{iopt} = [space.ncomp, rdim, npts_out];
        grad = true;
        
      case 'curl'
        if (ndim == 2 && rdim == 2)
          eu{iopt} = NaN (1, npts_out);
          eunum{iopt} = {1};
          eusize{iopt} = npts_out;
        elseif (ndim == 3 && rdim == 3)
          eu{iopt} = NaN (space.ncomp, rdim, npts_out);
          eunum{iopt} = {1:space.ncomp};
          eusize{iopt} = [rdim, npts_out];
        end
        curl = true;
        
      case 'divergence'
        eu{iopt} = NaN (1, npts_out);
        eunum{iopt} = {1};
        eusize{iopt} = npts_out;
        divergence = true;
    end
  end

  for ipt = 1:npts_on_F
    for idim = 1:ndim
      brk{idim} = zeta{idim}(ind{idim}(ipt):ind{idim}(ipt)+1);
    end
    msh = msh_cartesian (brk, num2cell(pts_param(:,ipt)), [], geometry, 'boundary', false);
    sp  = space.constructor (msh);

    msh_col = msh_evaluate_element_list (msh, 1);
    sp_col  = sp_evaluate_element_list (sp, msh_col, 'value', value, 'gradient', grad, ...
            'curl', curl, 'divergence', divergence);
    eu_aux = sp_eval_msh (u, sp_col, msh_col, options);

    for iopt = 1:nopts
      eu{iopt}(eunum{iopt}{:},inds_on_F(ipt)) = eu_aux{iopt};
    end  
  end
  for iopt = 1:nopts
    eu{iopt} = reshape (eu{iopt}, [eusize{iopt}, 1]); % The extra 1 makes things work also in 1D
  end

  if (nopts == 1)
    eu = eu{1};
  end

end
