% SP_EVAL_MSH: Evaluate a function, given by its degrees of freedom, at the points given by a msh object.
%
%   [eu, F] = sp_eval_msh (u, space, msh, [options], [lambda_lame, mu_lame]);
%
% INPUT:
%     
%     u:           vector of dof weights
%     space:       struct defining the discrete space (see for instance sp_scalar/sp_evaluate_col)
%     msh:         struct defining the points where to evaluate (see msh_cartesian/msh_evaluate_col)
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

  output_cell = true;
  if (nargin < 4)
    options = {'value'};
    output_cell = false;
  elseif (~iscell (options))
    options = {options};
    output_cell = false;
  end

  nopts = numel (options);
  eu       = cell (nopts, 1);
  wsize    = cell (nopts, 1);
  shp_size = cell (nopts, 1);
  eusize   = cell (nopts, 1);
  field    = cell (nopts, 1);

  for iopt = 1:nopts
    switch (lower (options{iopt}))
      case 'value'
        eu{iopt} = zeros (space.ncomp, msh.nqn, msh.nel);
        wsize{iopt} = [1 1];
        shp_size{iopt} = space.ncomp;
        field{iopt} = 'shape_functions';
        ind(iopt) = 3;
        eusize{iopt} = {1:space.ncomp, 1:msh.nqn};

      case 'gradient'
        eu{iopt} = zeros (space.ncomp, msh.rdim, msh.nqn, msh.nel);
        wsize{iopt} = [1 1 1];
        shp_size{iopt} = [space.ncomp, msh.rdim];
        field{iopt} = 'shape_function_gradients';
        ind(iopt) = 4;
        eusize{iopt} = {1:space.ncomp, 1:msh.rdim, 1:msh.nqn};
        
      case 'laplacian'
        eu{iopt} = zeros (space.ncomp, msh.nqn, msh.nel);
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
        field{iopt} = 'shape_function_curls';
        
      case 'divergence'
        if (space.ncomp == 1)
          error ('sp_eval_msh: the divergence is not computed for scalars')
        end
        eu{iopt} = zeros (msh.nqn, msh.nel);
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
        wsize{iopt} = {[1 1 1], 1};
        shp_size{iopt} = {[space.ncomp, msh.rdim], []};
        field{iopt} = {'shape_function_gradients', 'shape_function_divs'};
        indaux = [4, 2];
        eusize{iopt} = {1:space.ncomp, 1:msh.rdim, 1:msh.nqn};

      otherwise
        error ('sp_eval_msh: unknown option: %s', options{iopt})
    end
  end
  
  F = msh.geo_map;

  uc = zeros (size (space.connectivity));
  uc(space.connectivity~=0) = ...
        u(space.connectivity(space.connectivity~=0));
     
  for iopt = 1:nopts
    if (~strcmpi (options{iopt}, 'stress'))
      weight = reshape (uc, [wsize{iopt}, space.nsh_max, msh.nel]);
      space.(field{iopt}) = reshape (space.(field{iopt}), ...
                   [shp_size{iopt}, msh.nqn, space.nsh_max, msh.nel]);

      eu{iopt}(eusize{iopt}{:},:) = reshape (sum (bsxfun (@times, weight, ...
           space.(field{iopt})), ind(iopt)), [shp_size{iopt}, msh.nqn, msh.nel]);
    else
      for idim = 1:msh.rdim
        x{idim} = reshape (F(idim,:,:), msh.nqn, msh.nel);
      end
      mu_values = reshape (mu_lame (x{:}), 1, 1, msh.nqn, msh.nel);
      lambda_values = lambda_lame (x{:});
      for ii = 1:numel (wsize{iopt})
        weight = reshape (uc, [wsize{iopt}{ii}, space.nsh_max, msh.nel]);
        space.(field{iopt}{ii}) = reshape (space.(field{iopt}{ii}), ...
                   [shp_size{iopt}{ii}, msh.nqn, space.nsh_max, msh.nel]);
        aux{ii} = reshape (sum (bsxfun (@times, weight, ...
         space.(field{iopt}{ii})), indaux(ii)), [shp_size{iopt}{ii}, msh.nqn, msh.nel]);
      end
      eu{iopt}(eusize{iopt}{:}, :) = ...
          bsxfun (@times, mu_values, (aux{1} + permute (aux{1}, [2 1 3 4])));
      for idim = 1:msh.rdim
        eu{iopt}(idim,idim,:,:) = eu{iopt}(idim,idim,:,:) + ...
          reshape (lambda_values .* aux{2}, 1, 1, msh.nqn, msh.nel);
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
  
  if (~output_cell)
    eu = eu{1};
  end
  
end