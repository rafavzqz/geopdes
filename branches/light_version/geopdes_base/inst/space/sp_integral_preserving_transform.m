function sp = sp_integral_preserving_transform (sp, msh, value)

  if (value)
    jacdet = reshape (geopdes_det__ (msh.geo_map_jac), msh.nqn, 1, msh.nel);
    sp.shape_functions = bsxfun (@rdivide, sp.shape_functions, jacdet);
  end
  
  if (isfield (sp, 'shape_function_gradients'))
    rmfield (sp, 'shape_function_gradients')
  end
  if (isfield (sp, 'shape_function_hessians'))
    rmfield (sp, 'shape_function_hessians')
  end

end