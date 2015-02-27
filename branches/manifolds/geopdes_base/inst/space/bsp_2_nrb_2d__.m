% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
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

function sp = bsp_2_nrb_2d__ (sp, msh, W)

  W = repmat (reshape (W(sp.connectivity), 1, sp.nsh_max, msh.nel), [msh.nqn, 1, 1]);
  shape_functions = W.* sp.shape_functions;
  D = repmat (reshape (sum (shape_functions, 2), msh.nqn, 1, msh.nel), [1, sp.nsh_max, 1]);
  sp.shape_functions = shape_functions ./ D;

  if (isfield (sp, 'shape_function_gradients'))
    Bu = W .* reshape (sp.shape_function_gradients(1,:,:,:), [msh.nqn, sp.nsh_max, msh.nel]);
    Bv = W .* reshape (sp.shape_function_gradients(2,:,:,:), [msh.nqn, sp.nsh_max, msh.nel]);

    Du = repmat (reshape (sum (Bu, 2), msh.nqn, 1, msh.nel), [1, sp.nsh_max, 1]);
    Dv = repmat (reshape (sum (Bv, 2), msh.nqn, 1, msh.nel), [1, sp.nsh_max, 1]);

    Nu = (Bu - sp.shape_functions .* Du)./D;
    Nv = (Bv - sp.shape_functions .* Dv)./D;

    sp.shape_function_gradients (1,:,:,:) = Nu;
    sp.shape_function_gradients (2,:,:,:) = Nv;
  end
end
