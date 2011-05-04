% DO_SP_SCALAR_TO_VECTOR_3D__: * INTERNAL UNDOCUMENTED FUNCTION *
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

function sp = do_sp_scalar_to_vector_3d__ (spx, ndofx, spy, ndofy, spz, ndofz, msh, varargin)

gradient = true; divergence = false; curl = false;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_scalar_to_vector_3d: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'divergence'))
      divergence = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'curl'))
      curl = varargin {ii+1};
    else
      error ('sp_scalar_to_vector_3d: unknown option %s', varargin {ii});
    end
  end
end

if ((divergence || curl || gradient) ...
    && (~isfield (spx, 'shape_function_gradients') ...
        || (~isfield (spy, 'shape_function_gradients')) ...
        || (~isfield (spz, 'shape_function_gradients'))))
  error ('sp_scalar_to_vector_3d: component fields are missing required derivatives');
end

if ((~gradient) && (~divergence) && (~curl))
  if (isfield (spx, 'shape_function_gradients'))
    spx = rmfield (spx, 'shape_function_gradients');
  end
  if (isfield (spy, 'shape_function_gradients'))
    spy = rmfield (spy, 'shape_function_gradients');
  end
  if (isfield (spz, 'shape_function_gradients'))
    spz = rmfield (spz, 'shape_function_gradients');
  end  
end

sp.nsh_max      = spx.nsh_max + spy.nsh_max + spz.nsh_max;
sp.nsh          = spx.nsh + spy.nsh + spz.nsh;
sp.ndof         = spx.ndof + spy.ndof + spz.ndof;
sp.comp_dofs{1} = 1:ndofx;
sp.comp_dofs{2} = ndofx+(1:ndofy);
sp.comp_dofs{3} = ndofx+ndofy+(1:ndofz);
sp.connectivity = [spx.connectivity; spy.connectivity+spx.ndof; spz.connectivity+spx.ndof+spy.ndof];
sp.ncomp        = 3;
sp.ndof_dir     = [spx.ndof_dir; spy.ndof_dir; spz.ndof_dir];

sp.shape_functions = zeros (3, msh.nqn, sp.nsh_max, msh.nel);
sp.shape_functions(1,:, 1:spx.nsh_max,                         :)  = spx.shape_functions;
sp.shape_functions(2,:, spx.nsh_max+1:spx.nsh_max+spy.nsh_max, :)  = spy.shape_functions;
sp.shape_functions(3,:, spx.nsh_max+spy.nsh_max+1:sp.nsh_max,  :)  = spz.shape_functions;

if (isfield (spx,'dofs') && isfield (spy,'dofs') && isfield (spz,'dofs'))
  sp.dofs = [spx.dofs, spy.dofs+ndofx, spz.dofs+ndofx+ndofy];
  sp.comp_dofs{1} = intersect (sp.comp_dofs{1}, sp.dofs);
  sp.comp_dofs{2} = intersect (sp.comp_dofs{2}, sp.dofs);
  sp.comp_dofs{3} = intersect (sp.comp_dofs{3}, sp.dofs);
end

if (gradient)
  sp.shape_function_gradients = zeros (3, 3, msh.nqn, sp.nsh_max, msh.nel);
  sp.shape_function_gradients(1,:,:, 1:spx.nsh_max,                        :) = spx.shape_function_gradients;
  sp.shape_function_gradients(2,:,:, spx.nsh_max+1:spx.nsh_max+spy.nsh_max,:) = spy.shape_function_gradients;
  sp.shape_function_gradients(3,:,:, spx.nsh_max+spy.nsh_max+1:sp.nsh_max, :) = spz.shape_function_gradients;
end

if (divergence)
  sp.shape_function_divs = zeros (msh.nqn, sp.nsh_max, msh.nel);
  sp.shape_function_divs(:, 1:spx.nsh_max, :)                         = squeeze (spx.shape_function_gradients(1,:,:,:));
  sp.shape_function_divs(:, spx.nsh_max+1:spx.nsh_max+spy.nsh_max, :) = squeeze (spy.shape_function_gradients(2,:,:,:));
  sp.shape_function_divs(:, spx.nsh_max+spy.nsh_max+1:sp.nsh_max, :)  = squeeze (spz.shape_function_gradients(3,:,:,:));
end

if (curl)
  sp.shape_function_curls = zeros (3, msh.nqn, sp.nsh_max, msh.nel);
  sp.shape_function_curls(2,:,1:spx.nsh_max,:) = spx.shape_function_gradients(3,:,:,:);
  sp.shape_function_curls(3,:,1:spx.nsh_max,:) = -spx.shape_function_gradients(2,:,:,:);
  sp.shape_function_curls(1,:,spx.nsh_max+[1:spy.nsh_max],:) = -spy.shape_function_gradients(3,:,:,:);
  sp.shape_function_curls(3,:,spx.nsh_max+[1:spy.nsh_max],:) = spy.shape_function_gradients(1,:,:,:);
  sp.shape_function_curls(1,:,(spx.nsh_max+spy.nsh_max+1):sp.nsh_max,:) = spz.shape_function_gradients(2,:,:,:);
  sp.shape_function_curls(2,:,(spx.nsh_max+spy.nsh_max+1):sp.nsh_max,:) = -spz.shape_function_gradients(1,:,:,:);
end

sp.spfun = @(MSH) do_sp_scalar_to_vector_3d__ (spx.spfun (MSH), ndofx, spy.spfun (MSH), ndofy, spz.spfun (MSH), ndofz, MSH);
end
