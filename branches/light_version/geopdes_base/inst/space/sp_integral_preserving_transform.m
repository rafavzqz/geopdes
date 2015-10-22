function sp = sp_integral_preserving_transform (sp, msh, value)

  if (value)
    jacdet = reshape (geopdes_det__ (msh.geo_map_jac), msh.nqn, 1, msh.nel);
    sp.shape_functions = bsxfun (@rdivide, sp.shape_functions, jacdet);
    sp.shape_functions = sp.shape_functions;
    if (isfield (msh, 'side_number'))
      sp.shape_functions = sp.shape_functions * (-1)^msh.side_number;
    end
  end
  
  if (isfield (sp, 'shape_function_gradients'))
    sp = rmfield (sp, 'shape_function_gradients');
  end
  if (isfield (sp, 'shape_function_hessians'))
    sp = rmfield (sp, 'shape_function_hessians');
  end

end