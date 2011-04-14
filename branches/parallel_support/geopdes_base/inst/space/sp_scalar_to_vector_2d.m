% SP_SCALAR_TO_VECTOR_2D: Construct a 2D vector valued function space given the function space for each component.
%
%     spv = sp_scalar_to_vector_2d (spx, spy, msh, 'option1', value1, ...)
%
% INPUTS:
%
%     spx,spy:        function spaces for the x and y components
%     msh:            structure containing the domain partition and the quadrature rule (see msh_push_forward_2d)
%    'option', value: additional optional parameters, currently available options are:
%            
%              Name     |   Default value |  Meaning
%           ------------+-----------------+-----------
%            gradient   |      true       |  compute shape_function_gradients
%           divergence  |      false      |  compute shape_function_divs
%              curl     |      flase      |  compute shape_function_curls
%
% OUTPUT:
%
%    spv: structure representing the vector valued function space
%
%
% Copyright (C) 2010 Carlo de Falco
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.
 
% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.
% Author: Carlo de Falco <cdf AT users.sourceforge.net>
% Created: 2010-07-21

function sp = sp_scalar_to_vector_2d (spx, spy, msh, varargin)

sp = do_sp_scalar_to_vector__ (spx, spx.ndof, spy, spy.ndof, msh, varargin{:});

if (isfield (msh, 'boundary'))
  for iside = 1:numel(msh.boundary)
    sp.boundary(iside) = do_sp_scalar_to_vector__ (spx.boundary(iside), spx.ndof, ...
                                                   spy.boundary(iside), spy.ndof, ...
                                                   msh.boundary(iside), 'gradient', false);
  end
end

end

function sp = do_sp_scalar_to_vector__ (spx, ndofx, spy, ndofy, msh, varargin)

gradient = true; divergence = false; curl = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_scalar_to_vector_2d: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'divergence'))
      divergence = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'curl'))
      curl = varargin {ii+1};
    else
      error ('sp_scalar_to_vector_2d: unknown option %s', varargin {ii});
    end
  end
end

if ((divergence || curl || gradient) ... 
    && (~isfield (spx, 'shape_function_gradients') ...
        || (~isfield (spy, 'shape_function_gradients'))))
  error ('sp_scalar_to_vector_2d: component fields are missing required derivatives');
end

  
if ((~gradient) && (~divergence) && (~curl))
  if (isfield (spx, 'shape_function_gradients'))
    spx = rmfield (spx, 'shape_function_gradients');
  end
  if (isfield (spy, 'shape_function_gradients'))
    spy = rmfield (spy, 'shape_function_gradients');
  end
end

sp.nsh_max      = spx.nsh_max + spy.nsh_max;
sp.nsh          = spx.nsh + spy.nsh;
sp.ndof         = spx.ndof + spy.ndof;
sp.comp_dofs{1} = 1:ndofx;
sp.comp_dofs{2} = ndofx+(1:ndofy);
sp.connectivity = [spx.connectivity; spy.connectivity+spx.ndof];
sp.ncomp = 2;


sp.shape_functions = zeros (2, msh.nqn, sp.nsh_max, msh.nel);
sp.shape_functions(1,:,1:spx.nsh_max,:)              = spx.shape_functions;
sp.shape_functions(2,:, spx.nsh_max+1:sp.nsh_max,:)  = spy.shape_functions;

if (isfield (spx,'dofs') && isfield (spy,'dofs'))
  sp.dofs = [spx.dofs, spy.dofs+ndofx];
  sp.comp_dofs{1} = intersect (sp.comp_dofs{1}, sp.dofs);
  sp.comp_dofs{2} = intersect (sp.comp_dofs{2}, sp.dofs);
end

if (gradient)

  sp.shape_function_gradients = zeros (2, 2, msh.nqn, sp.nsh_max, msh.nel);
  sp.shape_function_gradients(1,:,:,1:spx.nsh_max,:)             = spx.shape_function_gradients;
  sp.shape_function_gradients(2,:,:, spx.nsh_max+1:sp.nsh_max,:) = spy.shape_function_gradients;

end

if (divergence)
  sp.shape_function_divs = zeros (msh.nqn, sp.nsh_max, msh.nel);
  sp.shape_function_divs(:, 1:spx.nsh_max, :)            = squeeze (spx.shape_function_gradients(1,:,:,:));
  sp.shape_function_divs(:, spx.nsh_max+1:sp.nsh_max, :) = squeeze (spy.shape_function_gradients(2,:,:,:));
end

if (curl)
  sp.shape_function_curls = zeros (msh.nqn, sp.nsh_max, msh.nel);
  sp.shape_function_curls(:, 1:spx.nsh_max, :)            = squeeze (-spx.shape_function_gradients(2,:,:,:));
  sp.shape_function_curls(:, spx.nsh_max+1:sp.nsh_max, :) = squeeze (spy.shape_function_gradients(1,:,:,:));
end

sp.spfun = @(MSH) do_sp_scalar_to_vector__ (spx.spfun (MSH), ndofx, spy.spfun (MSH), ndofy, MSH);
end
