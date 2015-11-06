function sp = sp_vector_grad_preserving_transform (sp, msh, value, gradient, curl, divergence)

  if (nargin < 3)
    value = true;
  end
  if (nargin < 4)
    gradient = false;
  end
  if (nargin < 5)
    curl = false;
  end
  if (nargin < 6)
    divergence = false;
  end

  if (gradient || curl || divergence)
    JinvT = geopdes_invT__ (msh.geo_map_jac);
    JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh.nqn, msh.nel]);
    for icomp = 1:sp.ncomp
      sp.shape_function_gradients(icomp,:,:,:,:) = geopdes_prod__ (JinvT, ...
          reshape (sp.shape_function_gradients(icomp,:,:,:,:), msh.rdim, msh.nqn, sp.nsh_max, msh.nel));
    end
    
    if (divergence)
      if (sp.ncomp == msh.rdim)
        sp.shape_function_divs = zeros (msh.nqn, sp.nsh_max, msh.nel);
        for icomp = 1:sp.ncomp
          sp.shape_function_divs = sp.shape_function_divs + ...
            reshape (sp.shape_function_gradients(icomp,icomp,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
        end
      else
        error('The divergence is not implemented in the case that rdim != ncomp')
      end
    elseif (isfield (sp, 'shape_function_divs'))
      sp = rmfield (sp, 'shape_function_divs');
    end

    if (curl)
      if (sp.ncomp == 2 && msh.rdim == 2)
        sp.shape_function_curls = reshape (sp.shape_function_gradients(2,1,:,:,:) - ...
			 	       sp.shape_function_gradients(1,2,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
      elseif (sp.ncomp == 3 && msh.rdim == 3)
        sp.shape_function_curls = zeros (sp.ncomp, msh.nqn, sp.nsh_max, msh.nel);
        for icomp = 1:sp.ncomp
          ind1 = mod(icomp,3) + 1;
          ind2 = mod(ind1, 3) + 1;
          sp.shape_function_curls(icomp,:,:,:) = reshape (sp.shape_function_gradients(ind2,ind1,:,:,:) - ...
                     sp.shape_function_gradients(ind1,ind2,:,:,:), 1, msh.nqn, sp.nsh_max, msh.nel);
        end
      else
        error('The curl is not implemented in the case that rdim != ncomp')
      end
    elseif (isfield (sp, 'shape_function_curls'))
      sp = rmfield (sp, 'shape_function_curls');
    end
    
    if (~gradient)
      sp = rmfield (sp, 'shape_function_gradients');
    end
  end

  if (~value && isfield (sp, 'shape_functions'))
    sp = rmfield (sp, 'shape_functions');
  end
  
end