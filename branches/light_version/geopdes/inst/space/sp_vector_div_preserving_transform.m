function sp = sp_vector_div_preserving_transform (sp, msh, value, gradient, curl, divergence)

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

  jacdet = reshape (geopdes_det__ (msh.geo_map_jac), msh.nqn, msh.nel);
  
  if (divergence  && isfield (sp, 'shape_function_divs'))
    jacdet = reshape (jacdet, msh.nqn, 1, msh.nel);
    sp.shape_function_divs = bsxfun (@rdivide, sp.shape_function_divs, jacdet);
  end

  if (gradient || curl || (divergence && ~isfield (sp, 'shape_function_divs')))
    shape_fun_grads = piola_transform_grad__ (sp, msh, sp.shape_functions, sp.shape_function_gradients);
  
    if (gradient)
      sp.shape_function_gradients = shape_fun_grads;
    end
    
    if (divergence && ~isfield (sp, 'shape_function_divs'))
      sp.shape_function_divs = zeros (msh.nqn, sp.nsh_max, msh.nel);
      for icomp = 1:sp.ncomp
        sp.shape_function_divs = sp.shape_function_divs + sp.shape_function_gradients(icomp,icomp,:,:,:);
      end
    end

    if (curl)
      if (msh.ndim ~= msh.rdim)
        error ('The computation of the curl has not been implemented for manifolds')
      end
      if (sp.ncomp == 2)
        sp.shape_function_curls = reshape (shape_fun_grads(2,1,:,:,:) - ...
	 	       shape_fun_grads(1,2,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
      elseif (sp.ncomp == 3)
        shape_fun_curls = zeros (sp.ncomp, msh.nqn, sp.nsh_max, msh.nel);
        for icomp = 1:sp.ncomp
          ind1 = mod(icomp,3) + 1;
          ind2 = mod(ind1, 3) + 1;
          shape_fun_curls(icomp,:,:,:) = reshape (shape_fun_grads(ind2,ind1,:,:,:) - ...
                     shape_fun_grads(ind1,ind2,:,:,:), 1, msh.nqn, sp.nsh_max, msh.nel);
        end
        sp.shape_function_curls = shape_fun_curls;
      end
    end
  end

% The transformation is applied to the value only at the end, to avoid
%  conflicts with the computation of the gradient

% This is still not working fine. One is forced to compute the boundary with
% msh_eval_boundary_side, instead of using the _tp operators.
  if (value)
    sp.shape_functions = geopdes_prod__ (msh.geo_map_jac, sp.shape_functions);
    jacdet = reshape (jacdet, 1, msh.nqn, 1, msh.nel);
    sp.shape_functions = bsxfun (@rdivide, sp.shape_functions, jacdet);
  end
  
  if (~gradient && isfield (sp, 'shape_function_gradients'))
    sp = rmfield (sp, 'shape_function_gradients');
  end

  if (~value && isfield (sp, 'shape_functions'))
    sp = rmfield (sp, 'shape_functions');
  end
  
end