% -*- INTERNAL UNDOCUMENTED FUNCTION -*-
%
% Copyright (C) 2010 Carlo de Falco
% Copyright (C) 2015 Rafael Vazquez
% Copyright (C) 2023 Pablo Antolin, Luca Coradello
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

function sp = bsp_2_nrb__ (sp, msh, W)

% sp is a struct (called from sp_evaluate_col)
  if (isstruct (sp))
    gradient = isfield (sp, 'shape_function_gradients');
    hessian  = isfield (sp, 'shape_function_hessians');
    third_derivative  = isfield (sp, 'shape_function_third_derivatives');
    fourth_derivative  = isfield (sp, 'shape_function_fourth_derivatives');
  end

  W = repmat (reshape (W(sp.connectivity), 1, sp.nsh_max, msh.nel), [msh.nqn, 1, 1]);
  shape_functions = W .* sp.shape_functions;
  D = repmat (reshape (sum (shape_functions, 2), msh.nqn, 1, msh.nel), [1, sp.nsh_max, 1]);
  sp.shape_functions = shape_functions ./ D;

  if (gradient)
    for idim = 1:msh.ndim
      Bdir = W .* reshape (sp.shape_function_gradients(idim,:,:,:), [msh.nqn, sp.nsh_max, msh.nel]);
      Ddir{idim} = repmat (reshape (sum (Bdir, 2), msh.nqn, 1, msh.nel), [1, sp.nsh_max, 1]);

      aux_grad{idim} = (Bdir - sp.shape_functions .* Ddir{idim}) ./ D;
      sp.shape_function_gradients(idim,:,:,:) = aux_grad{idim};
    end
  end

% With aux_grad I can avoid to reshape
  if (hessian)
    for idim = 1:msh.ndim
      for jdim = 1:msh.ndim
        Bdir2 = W .* reshape (sp.shape_function_hessians(idim,jdim,:,:,:), [msh.nqn, sp.nsh_max, msh.nel]);
        Ddir2 = repmat (reshape (sum (Bdir2, 2), msh.nqn, 1, msh.nel), [1, sp.nsh_max, 1]);
        sp.shape_function_hessians(idim,jdim,:,:,:) = (Bdir2 - sp.shape_functions .* Ddir2 ...
            - aux_grad{idim} .* Ddir{jdim} - aux_grad{jdim} .* Ddir{idim}) ./ D;
      end
    end
  end

  if (third_derivative)
    for idim = 1:msh.ndim
      for jdim = 1:msh.ndim
          for kdim = 1:msh.ndim
            Bdir3 = W .* reshape (sp.shape_function_third_derivatives(idim,jdim,kdim,:,:,:), [msh.nqn, sp.nsh_max, msh.nel]);
            Ddir3{idim,jdim,kdim} = repmat (reshape (sum (Bdir3, 2), msh.nqn, 1, msh.nel), [1, sp.nsh_max, 1]);
            aux_third{idim,jdim,kdim} = (Bdir3 - sp.shape_functions .* Ddir3{idim,jdim,kdim} ...
                - aux_grad{idim} .* Ddir2{jdim,kdim} - aux_grad{jdim} .* Ddir2{idim,kdim} - aux_grad{kdim} .* Ddir2{idim,jdim} ...
                - aux_hes{idim,jdim} .* Ddir{kdim} - aux_hes{idim,kdim} .* Ddir{jdim} - aux_hes{jdim,kdim} .* Ddir{idim}) ./ D;
            sp.shape_function_third_derivatives(idim,jdim,kdim,:,:,:) = aux_third{idim,jdim,kdim};
          end
      end
    end
  end
 
  if (fourth_derivative)
    for idim = 1:msh.ndim
      for jdim = 1:msh.ndim
          for kdim = 1:msh.ndim
              for ldim = 1:msh.ndim
                Bdir4 = W .* reshape (sp.shape_function_fourth_derivatives(idim,jdim,kdim,ldim,:,:,:), [msh.nqn, sp.nsh_max, msh.nel]);
                Ddir4 = repmat (reshape (sum (Bdir4, 2), msh.nqn, 1, msh.nel), [1, sp.nsh_max, 1]);
                sp.shape_function_fourth_derivatives(idim,jdim,kdim,ldim,:,:,:) = (Bdir4 - sp.shape_functions .* Ddir4 ...
                    - aux_grad{idim} .* Ddir3{jdim,kdim,ldim} - aux_grad{jdim} .* Ddir3{idim,kdim,ldim} - aux_grad{kdim} .* Ddir3{idim,jdim,ldim} ...
                    - aux_grad{ldim} .* Ddir3{idim,jdim,kdim} - aux_hes{idim,jdim} .* Ddir2{kdim,ldim} - aux_hes{idim,kdim} .* Ddir2{jdim,ldim} ...
                    - aux_hes{idim,ldim} .* Ddir2{jdim,kdim} - aux_hes{jdim,kdim} .* Ddir2{idim,ldim} - aux_hes{jdim,ldim} .* Ddir2{idim,kdim} ...
                    - aux_hes{kdim,ldim} .* Ddir2{idim,jdim} - aux_third{idim,jdim,kdim} .* Ddir{ldim} - aux_third{idim,jdim,ldim} .* Ddir{kdim} ...
                    - aux_third{idim,kdim,ldim} .* Ddir{jdim} - aux_third{jdim,kdim,ldim} .* Ddir{idim} ) ./ D;
              end
          end
      end
    end
  end

end
