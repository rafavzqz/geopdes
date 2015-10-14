function sp = sp_vector_curl_preserving_transform (sp, msh, value, curl)

  if (nargin < 3)
    value = true;
  end
  if (nargin < 4)
    curl = false;
  end

  [JinvT, jacdet] = geopdes_invT__ (msh.geo_map_jac);

  if (value)
    JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh.nqn, msh.nel]);
    sp.shape_functions = geopdes_prod__ (JinvT, sp.shape_functions);
  elseif (isfield (sp, 'shape_functions'))
    sp = rmfield (sp, 'shape_functions');
  end

  if (curl)
    if (sp.ncomp_param == 2)
      jacdet = reshape (jacdet, msh.nqn, 1, msh.nel);
      sp.shape_function_curls = bsxfun (@rdivide, sp.shape_function_curls, jacdet);
  
    elseif (sp.ncomp_param == 3)
      jacdet = reshape (jacdet, 1, msh.nqn, 1, msh.nel);
      sp.shape_function_curls = geopdes_prod__ (msh.geo_map_jac, sp.shape_function_curls);
      sp.shape_function_curls = bsxfun (@rdivide, sp.shape_function_curls, jacdet);
    end
  elseif (isfield (sp, 'shape_function_curls'))
    sp = rmfield (sp, 'shape_function_curls');
  end

  
  if (isfield (sp, 'shape_function_divs'))
    sp = rmfield (sp, 'shape_function_divs');
  end
    
  if (isfield (sp, 'shape_function_gradients'))
    sp = rmfield (sp, 'shape_function_gradients');
  end
  
end