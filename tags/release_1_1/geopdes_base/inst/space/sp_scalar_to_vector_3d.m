% SP_SCALAR_TO_VECTOR_3D: Construct a 3D vector valued function space given the function space for each component.
%
%     spv = sp_scalar_to_vector_3d (spx, spy, spz, msh, 'option1', value1, ...)
%
% INPUTS:
%
%     spx,spy,spz:    function spaces for the x, y and z components
%     msh:            structure containing the domain partition and the quadrature rule (see msh_push_forward_3d)
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
% Created: 2010-08-02

function sp = sp_scalar_to_vector_3d (spx, spy, spz, msh, varargin)

sp = do_sp_scalar_to_vector_3d__ (spx, spx.ndof, spy, spy.ndof, spz, spz.ndof, msh, varargin{:});

if (isfield (msh, 'boundary'))
  for iside = 1:numel(msh.boundary)
    sp.boundary(iside) = do_sp_scalar_to_vector_3d__ (spx.boundary(iside), spx.ndof, ...
                                                      spy.boundary(iside), spy.ndof, ...
                                                      spz.boundary(iside), spz.ndof, ...
                                                      msh.boundary(iside), 'gradient', false);
  end
end

end