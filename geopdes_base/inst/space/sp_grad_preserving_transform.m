function sp = sp_grad_preserving_transform (sp, msh, value, gradient, hessian, laplacian)

  if (hessian || laplacian)
    [Jinv, Jinv2] = geopdes_inv_der2__ (msh.geo_map_jac, msh.geo_map_der2);
    Jinv  = reshape (Jinv, [msh.ndim, msh.rdim, msh.nqn, 1, msh.nel]);
    JinvT = permute (Jinv, [2 1 3 4 5]);
    Jinv2 = reshape (Jinv2, [msh.ndim, msh.rdim, msh.rdim, msh.nqn, 1, msh.nel]);

    shg = reshape (sp.shape_function_gradients, [msh.ndim, 1, 1, msh.nqn, sp.nsh_max, msh.nel]);
    shh = reshape (sp.shape_function_hessians, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
    shape_function_hessians = zeros (msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);

    shh_size = [1, 1, msh.nqn, sp.nsh_max, msh.nel];
    for idim = 1:msh.rdim
      for jdim = 1:msh.rdim
        D2v_DF = sum (bsxfun(@times, shh, Jinv(:,idim,:,:,:)),1);
        DFT_D2v_DF = sum (bsxfun (@times, JinvT(jdim,:,:,:,:), D2v_DF), 2);
        Dv_D2F = sum (bsxfun (@times, shg, Jinv2(:,idim,jdim,:,:,:)), 1);

        shape_function_hessians(idim,jdim,:,:,:) = reshape (DFT_D2v_DF, shh_size) + ...
            reshape (Dv_D2F, shh_size);
      end
    end
      
    if (hessian)
      sp.shape_function_hessians = shape_function_hessians;
    end
      
    if (laplacian)
      sp.shape_function_laplacians = zeros (msh.nqn, sp.nsh_max, msh.nel);
      for idim = 1:msh.rdim
        sp.shape_function_laplacians = sp.shape_function_laplacians + ...
            reshape (shape_function_hessians(idim,idim,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
      end
    end
  end
  if (~hessian && isfield (sp, 'shape_function_hessians'))
      sp = rmfield (sp, 'shape_function_hessians');
  end
  if (~laplacian && isfield (sp, 'shape_function_laplacians'))
      sp = rmfield (sp, 'shape_function_laplacians');
  end

  if (gradient)
    JinvT = geopdes_invT__ (msh.geo_map_jac);
    JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh.nqn, msh.nel]);
    sp.shape_function_gradients = geopdes_prod__ (JinvT, sp.shape_function_gradients);
  elseif (isfield (sp, 'shape_function_gradients'))
    sp = rmfield (sp, 'shape_function_gradients');
  end

  if (~value && isfield (sp, 'shape_functions'))
    sp = rmfield (sp, 'shape_functions');
  end
  
end