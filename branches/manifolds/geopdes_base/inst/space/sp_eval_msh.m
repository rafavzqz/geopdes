% SP_EVAL_MSH: Evaluate a function, given by its degrees of freedom, at the points given by a msh object.
%
%   [eu, F] = sp_eval_msh (u, space, msh, [options], [lambda_lame, mu_lame]);
%
% INPUT:
%     
%     u:           vector of dof weights
%     space:       object defining the discrete space (see sp_bspline)
%     msh:         object defining the points where to evaluate (see msh_cartesian)
%     options:     cell array with the fields to plot
%                   accepted options are 'value' (default), 'gradient',
%                   and for vectors also 'curl', 'divergence', 'stress'
%     lambda_lame: function handle to the first Lame coefficient (only needed to compute 'stress')
%     mu_lame:     function handle for the second Lame coefficient (only needed to compute 'stress')
%
% OUTPUT:
%
%     eu: the function evaluated in the given points 
%     F:  grid points in the physical domain, that is, the mapped points
% 
% Copyright (C) 2009, 2010 Carlo de Falco
% Copyright (C) 2011, 2013, 2015 Rafael Vazquez
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

function [eu, F] = sp_eval_msh (u, space, msh, options, lambda_lame, mu_lame)

  if (nargin < 4)
    options = {'value'};
  elseif (~iscell (options))
    options = {options};
  end

  nopts = numel (options);
  F  = zeros (msh.rdim, msh.nqn, msh.nel);
  eu = cell (nopts, 1);
  wsize = cell (nopts, 1);

  value = false; grad = false; laplacian = false; curl = false; divergence = false;
  stress = false;
  
  for iopt = 1:nopts
    switch (lower (options{iopt}))
      case 'value'
        eu{iopt} = zeros (space.ncomp, msh.nqn, msh.nel);
        value = true;
        wsize{iopt} = [1 1];
        shp_size{iopt} = space.ncomp;
        field{iopt} = 'shape_functions';
        ind(iopt) = 3;
        eusize{iopt} = {1:space.ncomp, 1:msh.nqn};

      case 'gradient'
        eu{iopt} = zeros (space.ncomp, msh.rdim, msh.nqn, msh.nel);
        grad = true;
        wsize{iopt} = [1 1 1];
        shp_size{iopt} = [space.ncomp, msh.rdim];
        field{iopt} = 'shape_function_gradients';
        ind(iopt) = 4;
        eusize{iopt} = {1:space.ncomp, 1:msh.rdim, 1:msh.nqn};
        
      case 'laplacian'
        eu{iopt} = zeros (space.ncomp, msh.nqn, msh.nel);
        laplacian = true;
        wsize{iopt} = [1 1];
        shp_size{iopt} = space.ncomp;
        field{iopt} = 'shape_function_laplacians';
        ind(iopt) = 3;
        eusize{iopt} = {1:space.ncomp, 1:msh.nqn};
        
      case 'curl'
        if (space.ncomp == 1)
          error ('sp_eval_msh: the curl is not computed for scalars')
        end
        if (msh.ndim == 2 && msh.rdim == 2)
          eu{iopt} = zeros (msh.nqn, msh.nel);
          wsize{iopt} = 1;
          shp_size{iopt} = [];
          ind(iopt) = 2;
          eusize{iopt} = {1:msh.nqn};
        elseif (msh.ndim == 3 && msh.rdim == 3)
          eu{iopt} = zeros (msh.rdim, msh.nqn, msh.nel);
          wsize{iopt} = [1 1];
          shp_size{iopt} = space.ncomp;
          ind(iopt) = 3;
          eusize{iopt} = {1:msh.rdim, 1:msh.nqn};
        end
        curl = true;
        field{iopt} = 'shape_function_curls';
        
      case 'divergence'
        if (space.ncomp == 1)
          error ('sp_eval_msh: the divergence is not computed for scalars')
        end
        eu{iopt} = zeros (msh.nqn, msh.nel);
        divergence = true;
        wsize{iopt} = 1;
        shp_size{iopt} = [];
        field{iopt} = 'shape_function_divs';
        ind(iopt) = 2;
        eusize{iopt} = {1:msh.nqn};

      case 'stress'
        if (space.ncomp == 1)
          error ('sp_eval_msh: the stress is not computed for scalars')
        end
        if (nargin < 6)
          error ('sp_eval_msh: Lame coefficients missing')
        end
        eu{iopt} = zeros (space.ncomp, msh.rdim, msh.nqn, msh.nel);
        aux{1} = zeros (space.ncomp, msh.rdim, msh.nqn, msh.nel);
        aux{2} = zeros (msh.nqn, msh.nel);
        grad = true; divergence = true;
        wsize{iopt} = {[1 1 1], 1};
        shp_size{iopt} = {[space.ncomp, msh.rdim], []};
        field{iopt} = {'shape_function_gradients', 'shape_function_divs'};
        indaux = [4, 2];
        eusize{iopt} = {1:space.ncomp, 1:msh.rdim, 1:msh.nqn};

    end
  end
  
  

  for iel = 1:msh.nel_dir(1)
    msh_col = msh_evaluate_col (msh, iel);
    if (space.ncomp == 1)
      sp_col  = sp_evaluate_col (space, msh_col, 'value', value, 'gradient', grad, ...
          'laplacian', laplacian);
    else
      sp_col  = sp_evaluate_col (space, msh_col, 'value', value, 'gradient', grad, ...
          'curl', curl, 'divergence', divergence);
    end

    F(:,:,msh_col.elem_list) = msh_col.geo_map;

    uc_iel = zeros (size (sp_col.connectivity));
    uc_iel(sp_col.connectivity~=0) = ...
          u(sp_col.connectivity(sp_col.connectivity~=0));

    for iopt = 1:nopts
      if (~strcmpi (options{iopt}, 'stress'))
        weight = reshape (uc_iel, [wsize{iopt}, sp_col.nsh_max, msh_col.nel]);
        sp_col.(field{iopt}) = reshape (sp_col.(field{iopt}), ...
                     [shp_size{iopt}, msh_col.nqn, sp_col.nsh_max, msh_col.nel]);

        eu{iopt}(eusize{iopt}{:},msh_col.elem_list) = reshape (sum (bsxfun (@times, weight, ...
             sp_col.(field{iopt})), ind(iopt)), [shp_size{iopt}, msh_col.nqn, msh_col.nel]);
       
      else
        for idim = 1:msh.rdim
          x{idim} = reshape (F(idim,:,msh_col.elem_list), msh_col.nqn, msh_col.nel);
        end
        mu_col = reshape (mu_lame (x{:}), 1, 1, msh_col.nqn, msh_col.nel);
        lambda_col = lambda_lame (x{:});

        for ii = 1:numel (wsize{iopt})
          weight = reshape (uc_iel, [wsize{iopt}{ii}, sp_col.nsh_max, msh_col.nel]);
          sp_col.(field{iopt}{ii}) = reshape (sp_col.(field{iopt}{ii}), ...
                     [shp_size{iopt}{ii}, msh_col.nqn, sp_col.nsh_max, msh_col.nel]);
          aux{ii} = reshape (sum (bsxfun (@times, weight, ...
           sp_col.(field{iopt}{ii})), indaux(ii)), [shp_size{iopt}{ii}, msh_col.nqn, msh_col.nel]);
        end
        eu{iopt}(eusize{iopt}{:}, msh_col.elem_list) = ...
            bsxfun (@times, mu_col, (aux{1} + permute (aux{1}, [2 1 3 4])));
        for idim = 1:msh.rdim
          eu{iopt}(idim,idim,:,msh_col.elem_list) = eu{iopt}(idim,idim,:,msh_col.elem_list) + ...
            reshape (lambda_col .* aux{2}, 1, 1, msh_col.nqn, msh_col.nel);
        end
      end
       
    end
  end

  if (space.ncomp == 1)
    for iopt = 1:nopts
      switch (lower (options{iopt}))
        case {'value', 'laplacian'}
          eu{iopt} = reshape (eu{iopt}, msh.nqn, msh.nel);
        case {'gradient'}
          eu{iopt} = reshape (eu{iopt}, msh.rdim, msh.nqn, msh.nel);
      end
    end
  end
  
  % For compatibility with previous versions
  if (nopts == 1)
    eu = eu{1};
  end
  
end