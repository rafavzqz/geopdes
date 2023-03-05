% SP_GRAD_PRESERVING_TRANSFORM: apply the grad-preserving transform to the functions in the parametric domain
%
%     sp = sp_grad_preserving_transform (space, msh, [value, gradient, hessian, laplacian])
%
% INPUTS:
%     
%    space:   structure with the information in the parametric domain (see sp_scalar/sp_evaluate_col)
%    msh:     msh structure containing the information of the parametrization
%              in the points where basis functions have to be computed (see msh_cartesian/msh_evaluate_col)
%    value, gradient, hessian, laplacian: additional optional parameters, either true or false
%            
%              Name             |   Default value |  Meaning
%           --------------------+-----------------+----------------------------------
%            value              |      true       |  compute shape_functions
%            gradient           |      false      |  compute shape_function_gradients
%            hessian            |      false      |  compute shape_function_hessians
%            laplacian          |      false      |  compute shape_function_laplacians
%            third_derivative   |      false      |  compute shape_function_third_derivatives
%            fourth_derivative  |      false      |  compute shape_function_fourth_derivatives
%            bilaplacian        |      false      |  compute shape_function_bilaplacians
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%              (see the article for a detailed description)
%
%    FIELD_NAME      (SIZE)                                 DESCRIPTION
%    ncomp           (scalar)                               number of components of the functions of the space (actually, 1)
%    ndof            (scalar)                               total number of degrees of freedom
%    ndof_dir        (ncomp x ndim matrix)                  for each component, number of degrees of freedom along each direction
%    nsh_max         (scalar)                               maximum number of shape functions per element
%    nsh             (1 x msh.nel vector)                   actual number of shape functions per each element
%    connectivity    (nsh_max x msh.nel vector)             indices of basis functions that do not vanish in each element
%    shape_functions (msh.nqn x nsh_max x msh.nel)          basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%       (rdim x msh.nqn x nsh_max x msh.nel)                basis function gradients evaluated at each quadrature node in each element
%    shape_function_hessians
%       (rdim x rdim x msh.nqn x nsh_max x msh.nel)         basis function hessians evaluated at each quadrature node in each element
%    shape_function_laplacians 
%       (msh.nqn x nsh_max x msh.nel)                       basis functions laplacians evaluated at each quadrature node in each element
%    shape_function_third_derivatives
%       (rdim x rdim x rdim x msh.nqn x nsh_max x msh.nel)  basis function third derivatives evaluated at each quadrature node in each element
%    shape_function_fourth_derivatives
%       (rdim x rdim x rdim x rdim x msh.nqn x nsh_max x msh.nel) basis function fourth derivatives evaluated at each quadrature node in each element
%    shape_function_bilaplacians
%       (msh.nqn x nsh_max x msh.nel)               basis function bilaplacians evaluated at each quadrature node in each element
%
% Copyright (C) 2015 Rafael Vazquez
% Copyright (C) 2023 Pablo Antolin, Luca Coradello
%
%    This program is free software: you can redistribute it and/or modify
%    it under the terms of the GNU General Public License as published by
%    the Free Software Foundation, either version 3 of the License, or
%    (at your option) any later version.

%    This program is distributed in the hope that it will be useful,
%    but WITHOUT ANY WARRANTY; without even the implied warranty of
%    MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
%    GNU General Public License for more details.
%
%    You should have received a copy of the GNU General Public License
%    along with this program.  If not, see <http://www.gnu.org/licenses/>.

function sp = sp_grad_preserving_transform (sp, msh, value, gradient, hessian, laplacian, third_derivative, fourth_derivative, bilaplacian)

  if (nargin < 3 || isempty (value))
    value = true;
  end
  if (nargin < 4 || isempty (gradient))
    gradient = false;
  end
  if (nargin < 5 || isempty (hessian))
    hessian = false;
  end
  if (nargin < 6 || isempty (laplacian))
    laplacian = false;
  end
  if (nargin < 7 || isempty (third_derivative))
    third_derivative = false;
  end
  if (nargin < 8 || isempty (fourth_derivative))
    fourth_derivative = false;
  end
  if (nargin < 9 || isempty (bilaplacian))
    bilaplacian = false;
  end

  if (~value && isfield (sp, 'shape_functions'))
    sp = rmfield (sp, 'shape_functions');
  end
  
  if (gradient || laplacian || hessian || third_derivative || fourth_derivative || bilaplacian)
    JinvT = geopdes_invT__ (msh.geo_map_jac);
    JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh.nqn, msh.nel]);
    shape_gradients = geopdes_prod__ (JinvT, sp.shape_function_gradients);
    if(gradient)
        sp.shape_function_gradients = shape_gradients;
    end
  elseif (isfield (sp, 'shape_function_gradients'))
    sp = rmfield (sp, 'shape_function_gradients');
  end
  
  if (hessian || laplacian || third_derivative || fourth_derivative || bilaplacian)
    JinvT = geopdes_invT__ (msh.geo_map_jac);
    
    JinvT1 = reshape (JinvT, [msh.rdim, 1, msh.ndim, 1, msh.nqn, 1, msh.nel]);
    JinvT2 = reshape (JinvT, [1, msh.rdim, 1, msh.ndim, msh.nqn, 1, msh.nel]);

    J2 = reshape (msh.geo_map_der2, [msh.rdim, msh.ndim, msh.ndim, msh.nqn, 1, msh.nel]);
    
    shape_gradients = reshape (shape_gradients, [msh.rdim, 1, 1, msh.nqn, sp.nsh_max, msh.nel]);
        
    shh = reshape (sp.shape_function_hessians, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
    
    % See Appendix paper H. Casquero et al. 2016 to understand how to compute higher order derivatives 
    tmp_product = sum (bsxfun (@times, shape_gradients, J2), 1);

    tmp_product = reshape (tmp_product, [msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel]);

    tmp_subtraction = shh - tmp_product;

    tmp_subtraction = reshape (tmp_subtraction, [1, 1, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel]);

    shape_hessians = sum ( sum (bsxfun (@times, bsxfun(@times, JinvT1, tmp_subtraction), JinvT2 ), 3), 4);

    shape_hessians = reshape (shape_hessians, [msh.rdim, msh.rdim, msh.nqn, sp.nsh_max, msh.nel]);
  
    if (hessian)
      sp.shape_function_hessians = shape_hessians;
    end
      
    if (laplacian)
      sp.shape_function_laplacians = zeros (msh.nqn, sp.nsh_max, msh.nel);
      for idim = 1:msh.rdim
        sp.shape_function_laplacians = sp.shape_function_laplacians + ...
            reshape (shape_hessians(idim,idim,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
      end
    end
  end
  if (~hessian && isfield (sp, 'shape_function_hessians'))
      sp = rmfield (sp, 'shape_function_hessians');
  end
  if (~laplacian && isfield (sp, 'shape_function_laplacians'))
      sp = rmfield (sp, 'shape_function_laplacians');
  end

  if (third_derivative || fourth_derivative || bilaplacian)     % See Appendix paper H. Casquero et al. to understand how to compute higher order derivatives 

    JinvT = geopdes_invT__ (msh.geo_map_jac);
    
    JinvT1 = reshape (JinvT, [msh.rdim, 1, 1, msh.ndim, 1, 1, msh.nqn, 1, msh.nel]);
    JinvT2 = reshape (JinvT, [1, msh.rdim, 1, 1, msh.ndim, 1, msh.nqn, 1, msh.nel]);
    JinvT3 = reshape (JinvT, [1, 1, msh.rdim, 1, 1, msh.ndim, msh.nqn, 1, msh.nel]);
    
    
    J_1 = reshape(msh.geo_map_jac, [msh.rdim, 1, msh.ndim, 1, 1, msh.nqn, 1, msh.nel]);
    J_2 = reshape(msh.geo_map_jac, [1, msh.rdim, 1, msh.ndim, 1, msh.nqn, 1, msh.nel]);
    J_3 = reshape(msh.geo_map_jac, [1, msh.rdim, 1, 1, msh.ndim, msh.nqn, 1, msh.nel]);
   
    J2_1 = reshape (msh.geo_map_der2, [msh.rdim, 1, msh.ndim, 1, msh.ndim, msh.nqn, 1, msh.nel]);
    J2_2 = reshape (msh.geo_map_der2, [1, msh.rdim, 1, msh.ndim, msh.ndim, msh.nqn, 1, msh.nel]);
    J2_3 = reshape (msh.geo_map_der2, [msh.rdim, 1, msh.ndim, msh.ndim, 1, msh.nqn, 1, msh.nel]);
    
    shape_hessians = reshape(shape_hessians, [msh.rdim, msh.rdim, 1, 1, 1, msh.nqn, sp.nsh_max, msh.nel]);
    
    tmp_product_1 = sum ( sum (bsxfun (@times, bsxfun(@times, shape_hessians, J_2), J2_1 ), 2), 1);
   
    tmp_product_2 = sum ( sum (bsxfun (@times, bsxfun(@times, shape_hessians, J_1), J2_2 ), 1), 2);
    
    tmp_product_3 = sum ( sum (bsxfun (@times, bsxfun(@times, shape_hessians, J_3), J2_3 ), 2), 1);
    
    tmp_sum = reshape((tmp_product_1 + tmp_product_2 + tmp_product_3), [msh.ndim, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel]);
    
    J3 = reshape (msh.geo_map_der3, [msh.rdim, msh.ndim, msh.ndim, msh.ndim, msh.nqn, 1, msh.nel]);

    shape_gradients = reshape (shape_gradients, [msh.rdim, 1, 1, 1, msh.nqn, sp.nsh_max, msh.nel]);
     
    tmp_product_4 = sum (bsxfun (@times, shape_gradients, J3), 1);

    tmp_product_4 = reshape (tmp_product_4, [msh.ndim, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel]);
    
    shtd = reshape(sp.shape_function_third_derivatives, [msh.ndim, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel]);
    
    tmp_subtraction = shtd - tmp_product_4 - tmp_sum;

    tmp_subtraction = reshape (tmp_subtraction, [1, 1, 1, msh.ndim, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel]);

    shape_third_derivatives = sum ( sum ( sum ( bsxfun(@times, bsxfun(@times, bsxfun(@times, tmp_subtraction, JinvT1), JinvT2 ), JinvT3), 4), 5), 6);

    shape_third_derivatives = reshape (shape_third_derivatives, [msh.rdim, msh.rdim, msh.rdim, msh.nqn, sp.nsh_max, msh.nel]);
      
    if (third_derivative)
      sp.shape_function_third_derivatives = shape_third_derivatives;
    end      
  end
  if (~third_derivative && isfield (sp, 'shape_function_third_derivatives'))
      sp = rmfield (sp, 'shape_function_third_derivatives');
  end

  if (fourth_derivative || bilaplacian)     % See Appendix paper H. Casquero et al. to understand how to compute higher order derivatives 

    JinvT = geopdes_invT__ (msh.geo_map_jac);
    
    JinvT1 = reshape (JinvT, [msh.rdim, 1, 1, 1, msh.ndim, 1, 1, 1, msh.nqn, 1, msh.nel]);
    JinvT2 = reshape (JinvT, [1, msh.rdim, 1, 1, 1, msh.ndim, 1, 1, msh.nqn, 1, msh.nel]);
    JinvT3 = reshape (JinvT, [1, 1, msh.rdim, 1, 1, 1, msh.ndim, 1, msh.nqn, 1, msh.nel]);
    JinvT4 = reshape (JinvT, [1, 1, 1, msh.rdim, 1, 1, 1, msh.ndim, msh.nqn, 1, msh.nel]);
    
    
    J_1 = reshape(msh.geo_map_jac, [msh.rdim, 1, 1, msh.ndim, 1, 1, 1, msh.nqn, 1, msh.nel]); % m,alpha
    J_2 = reshape(msh.geo_map_jac, [1, 1, msh.rdim, 1, 1, msh.ndim, 1, msh.nqn, 1, msh.nel]); % q,gamma
    J_3 = reshape(msh.geo_map_jac, [1, msh.rdim, 1, 1, msh.ndim, 1, 1, msh.nqn, 1, msh.nel]); % n,beta
  
    J2_1 = reshape (msh.geo_map_der2, [msh.rdim, 1, 1, msh.ndim, 1, 1, msh.ndim, msh.nqn, 1, msh.nel]); % m,alpha delta
    J2_2 = reshape (msh.geo_map_der2, [1, msh.rdim, 1, 1, msh.ndim, 1, msh.ndim, msh.nqn, 1, msh.nel]); % n,beta delta
    J2_3 = reshape (msh.geo_map_der2, [1, 1, msh.rdim, 1, 1, msh.ndim, msh.ndim, msh.nqn, 1, msh.nel]); % q,gamma delta
    
    shape_third_derivatives = reshape(shape_third_derivatives, [msh.rdim, msh.rdim, msh.rdim, 1, 1, 1, 1, msh.nqn, sp.nsh_max, msh.nel]);
    
    tmp_product_1 = sum ( sum ( sum (bsxfun (@times, bsxfun (@times, bsxfun(@times, shape_third_derivatives, J_3), J_2 ), J2_1 ), 2), 3), 1);
   
    tmp_product_2 = sum ( sum ( sum (bsxfun (@times, bsxfun (@times, bsxfun(@times, shape_third_derivatives, J_1), J_2 ), J2_2 ), 1), 3), 2);
    
    tmp_product_3 = sum ( sum ( sum (bsxfun (@times, bsxfun (@times, bsxfun(@times, shape_third_derivatives, J_1), J_3 ), J2_3 ), 1), 2), 3);

    clear J_1 J_2 J_3 J2_1 J2_2 J2_3                

    J_1 = reshape(msh.geo_map_jac, [msh.rdim, 1, 1, msh.ndim, 1, 1, 1, msh.nqn, 1, msh.nel]); % m,alpha
    J_2 = reshape(msh.geo_map_jac, [1, msh.rdim, 1, 1, 1, msh.ndim, 1, msh.nqn, 1, msh.nel]); % n,gamma
    J_3 = reshape(msh.geo_map_jac, [1, msh.rdim, 1, 1, msh.ndim, 1, 1, msh.nqn, 1, msh.nel]); % n,beta
    J_4 = reshape(msh.geo_map_jac, [1, 1, msh.rdim, 1, 1, 1, msh.ndim, msh.nqn, 1, msh.nel]); % q,delta
    
    J2_1 = reshape (msh.geo_map_der2, [msh.rdim, 1, 1, msh.ndim, msh.ndim, 1, 1, msh.nqn, 1, msh.nel]); % m,alpha beta
    J2_2 = reshape (msh.geo_map_der2, [msh.rdim, 1, 1, msh.ndim, 1, msh.ndim, 1, msh.nqn, 1, msh.nel]); % m,alpha gamma
    J2_3 = reshape (msh.geo_map_der2, [1, msh.rdim, 1, 1, msh.ndim, msh.ndim, 1, msh.nqn, 1, msh.nel]); % n,beta gamma
      
    tmp_product_4 = sum ( sum ( sum (bsxfun (@times, bsxfun (@times, bsxfun(@times, shape_third_derivatives, J_4), J2_2 ), J_3 ), 3), 1), 2);
   
    tmp_product_5 = sum ( sum ( sum (bsxfun (@times, bsxfun (@times, bsxfun(@times, shape_third_derivatives, J_4), J2_3 ), J_1 ), 3), 2), 1);
    
    tmp_product_6 = sum ( sum ( sum (bsxfun (@times, bsxfun (@times, bsxfun(@times, shape_third_derivatives, J_4), J2_1 ), J_2 ), 3), 1), 2);
    
    tmp_sum_1 = reshape((tmp_product_1 + tmp_product_2 + tmp_product_3 + tmp_product_4 + tmp_product_5 + tmp_product_6 ),...
                        [msh.ndim, msh.ndim, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel]);
  
    clear tmp_product_1 tmp_product_2 tmp_product_3 tmp_product_4 tmp_product_5 tmp_product_6                
    clear J_1 J_2 J_3 J_4 J2_1 J2_2 J2_3                
                    
    J3_1 = reshape (msh.geo_map_der3, [msh.rdim, 1, msh.ndim, 1, msh.ndim, msh.ndim, msh.nqn, 1, msh.nel]); % m, alpha gamma delta
    J3_2 = reshape (msh.geo_map_der3, [1, msh.rdim, 1, msh.ndim, msh.ndim, msh.ndim, msh.nqn, 1, msh.nel]); % n, beta gamma delta
    J3_3 = reshape (msh.geo_map_der3, [msh.rdim, 1, msh.ndim, msh.ndim, msh.ndim, 1, msh.nqn, 1, msh.nel]); % m, alpha beta gamma 
    J3_4 = reshape (msh.geo_map_der3, [msh.rdim, 1, msh.ndim, msh.ndim, 1, msh.ndim, msh.nqn, 1, msh.nel]); % m, alpha beta delta

    J2_1 = reshape (msh.geo_map_der2, [msh.rdim, 1, msh.ndim, 1, msh.ndim, 1, msh.nqn, 1, msh.nel]); % m,alpha gamma
    J2_2 = reshape (msh.geo_map_der2, [1, msh.rdim, 1, msh.ndim, 1, msh.ndim, msh.nqn, 1, msh.nel]); % n,beta delta
    
    J2_3 = reshape (msh.geo_map_der2, [msh.rdim, 1, msh.ndim, 1, 1, msh.ndim, msh.nqn, 1, msh.nel]); % m,alpha delta
    J2_4 = reshape (msh.geo_map_der2, [1, msh.rdim, 1, msh.ndim, msh.ndim, 1, msh.nqn, 1, msh.nel]); % n,beta gamma

    J2_5 = reshape (msh.geo_map_der2, [msh.rdim, 1, msh.ndim, msh.ndim, 1, 1, msh.nqn, 1, msh.nel]); % m,alpha beta
    J2_6 = reshape (msh.geo_map_der2, [1, msh.rdim, 1, 1, msh.ndim, msh.ndim, msh.nqn, 1, msh.nel]); % n,gamma delta
 
    J_1 = reshape(msh.geo_map_jac, [1, msh.rdim, 1, msh.ndim, 1, 1, msh.nqn, 1, msh.nel]); % n,beta
    J_2 = reshape(msh.geo_map_jac, [msh.rdim, 1, msh.ndim, 1, 1, 1, msh.nqn, 1, msh.nel]); % m,alpha
    J_3 = reshape(msh.geo_map_jac, [1, msh.rdim, 1, 1, msh.ndim, 1, msh.nqn, 1, msh.nel]); % n,gamma
    J_4 = reshape(msh.geo_map_jac, [1, msh.rdim, 1, 1, 1, msh.ndim, msh.nqn, 1, msh.nel]); % n,delta
    
    shape_hessians = reshape(shape_hessians, [msh.rdim, msh.rdim, 1, 1, 1, 1, msh.nqn, sp.nsh_max, msh.nel]);

    tmp_product_1 = sum ( sum (bsxfun (@times, bsxfun(@times, shape_hessians, J3_1), J_1 ), 1), 2);
   
    tmp_product_2 = sum ( sum (bsxfun (@times, bsxfun(@times, shape_hessians, J3_2), J_2 ), 2), 1);
    
    tmp_product_3 = sum ( sum (bsxfun (@times, bsxfun(@times, shape_hessians, J3_3), J_4 ), 1), 2);

    tmp_product_4 = sum ( sum (bsxfun (@times, bsxfun(@times, shape_hessians, J3_4), J_3 ), 1), 2);
   
    tmp_product_5 = sum ( sum (bsxfun (@times, bsxfun(@times, shape_hessians, J2_1), J2_2 ), 1), 2);
    
    tmp_product_6 = sum ( sum (bsxfun (@times, bsxfun(@times, shape_hessians, J2_3), J2_4 ), 1), 2);

    tmp_product_7 = sum ( sum (bsxfun (@times, bsxfun(@times, shape_hessians, J2_5), J2_6 ), 1), 2);

    tmp_sum_2 = reshape((tmp_product_1 + tmp_product_2 + tmp_product_3 + tmp_product_4 + tmp_product_5 + tmp_product_6 + tmp_product_7),...
                        [msh.ndim, msh.ndim, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel]);
    
    clear tmp_product_1 tmp_product_2 tmp_product_3 tmp_product_4 tmp_product_5 tmp_product_6 tmp_product_7                
        
    shape_gradients = reshape (shape_gradients, [msh.rdim, 1, 1, 1, 1, msh.nqn, sp.nsh_max, msh.nel]);

    J4 = reshape (msh.geo_map_der4, [msh.rdim, msh.ndim, msh.ndim, msh.ndim, msh.ndim, msh.nqn, 1, msh.nel]);
     
    tmp_product = sum (bsxfun (@times, shape_gradients, J4), 1);

    tmp_product = reshape (tmp_product, [msh.ndim, msh.ndim, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel]);
    
    shfd = reshape(sp.shape_function_fourth_derivatives, [msh.ndim, msh.ndim, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel]);
    
    tmp_subtraction = shfd - tmp_sum_1 - tmp_sum_2 - tmp_product;

    tmp_subtraction = reshape (tmp_subtraction, [1, 1, 1, 1, msh.ndim, msh.ndim, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel]);

    shape_fourth_derivatives = sum ( sum ( sum ( sum ( bsxfun(@times, bsxfun(@times, bsxfun(@times, bsxfun(@times, ...
                                tmp_subtraction, JinvT1), JinvT2 ), JinvT3), JinvT4) , 5), 6), 7), 8);

    shape_fourth_derivatives = reshape (shape_fourth_derivatives, [msh.rdim, msh.rdim, msh.rdim, msh.rdim, msh.nqn, sp.nsh_max, msh.nel]);

    clear tmp_product tmp_subtraction           

    if (fourth_derivative)
      sp.shape_function_fourth_derivatives = shape_fourth_derivatives;
    end 
    if (bilaplacian)
      sp.shape_function_bilaplacians = zeros (msh.nqn, sp.nsh_max, msh.nel);
      for idim = 1:msh.rdim
        for jdim = 1:msh.rdim  
            sp.shape_function_bilaplacians = sp.shape_function_bilaplacians + ...
                reshape (shape_fourth_derivatives(idim,idim,jdim,jdim,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
        end
      end
    end    
  end
  if (~fourth_derivative && isfield (sp, 'shape_function_fourth_derivatives'))
      sp = rmfield (sp, 'shape_function_fourth_derivatives');
  end
  if (~bilaplacian && isfield (sp, 'shape_function_bilaplacians'))
      sp = rmfield (sp, 'shape_function_bilaplacians');
  end  
  
%   if (hessian || laplacian)
%     [Jinv, Jinv2] = geopdes_inv_der2__ (msh.geo_map_jac, msh.geo_map_der2);
%     Jinv  = reshape (Jinv, [msh.ndim, msh.rdim, msh.nqn, 1, msh.nel]);
%     JinvT = permute (Jinv, [2 1 3 4 5]);
%     Jinv2 = reshape (Jinv2, [msh.ndim, msh.rdim, msh.rdim, msh.nqn, 1, msh.nel]);
% 
%     shg = reshape (sp.shape_function_gradients, [msh.ndim, 1, 1, msh.nqn, sp.nsh_max, msh.nel]);
%     shh = reshape (sp.shape_function_hessians, msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
%     shape_function_hessians = zeros (msh.ndim, msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
% 
%     shh_size = [1, 1, msh.nqn, sp.nsh_max, msh.nel];
%     for idim = 1:msh.rdim
%       for jdim = 1:msh.rdim
%         D2v_DF = sum (bsxfun(@times, shh, Jinv(:,idim,:,:,:)),1);
%         DFT_D2v_DF = sum (bsxfun (@times, JinvT(jdim,:,:,:,:), D2v_DF), 2);
%         Dv_D2F = sum (bsxfun (@times, shg, Jinv2(:,idim,jdim,:,:,:)), 1);
% 
%         shape_function_hessians(idim,jdim,:,:,:) = reshape (DFT_D2v_DF, shh_size) + ...
%             reshape (Dv_D2F, shh_size);
%       end
%     end
%       
%     if (hessian)
%       sp.shape_function_hessians = shape_function_hessians;
%     end
%       
%     if (laplacian)
%       sp.shape_function_laplacians = zeros (msh.nqn, sp.nsh_max, msh.nel);
%       for idim = 1:msh.rdim
%         sp.shape_function_laplacians = sp.shape_function_laplacians + ...
%             reshape (shape_function_hessians(idim,idim,:,:,:), msh.nqn, sp.nsh_max, msh.nel);
%       end
%     end
%   end
%   if (~hessian && isfield (sp, 'shape_function_hessians'))
%       sp = rmfield (sp, 'shape_function_hessians');
%   end
%   if (~laplacian && isfield (sp, 'shape_function_laplacians'))
%       sp = rmfield (sp, 'shape_function_laplacians');
%   end
% 
%   if (gradient)
%     JinvT = geopdes_invT__ (msh.geo_map_jac);
%     JinvT = reshape (JinvT, [msh.rdim, msh.ndim, msh.nqn, msh.nel]);
%     sp.shape_function_gradients = geopdes_prod__ (JinvT, sp.shape_function_gradients);
%   elseif (isfield (sp, 'shape_function_gradients'))
%     sp = rmfield (sp, 'shape_function_gradients');
%   end
    
end

%!test
%! degree     = [5 5];       % Degree of the splines
%! regularity = [4 4];       % Regularity of the splines
%! nsub       = [3 3];       % Number of subdivisions
%! nquad      = [1 1];       % Points for the Gaussian quadrature rule
%! geo_name = 'geo_ring_1eighth.txt';
%! geometry = geo_load (geo_name);
%! degelev  = max (degree - (geometry.nurbs.order-1), 0);
%! nurbs    = nrbdegelev(geometry.nurbs, degelev);
%! [rknots, zeta, nknots] = kntrefine (nurbs.knots, nsub-1, nurbs.order-1, regularity);
%!
%! nurbs = nrbkntins (nurbs, nknots);
%! geometry = geo_load (nurbs);
%! rule     = msh_gauss_nodes (nquad);
%! [qn, qw] = msh_set_quad_nodes (zeta, rule);
%! msh      = msh_cartesian (zeta, qn, qw, geometry);
%! space  = sp_nurbs (geometry.nurbs, msh);
%!
%! msh_col = msh_evaluate_col (msh, 1);
%!
%! fourth_param = true;
%! third_param = true;
%! hessian_param = true;
%! grad_param = true;
%! value_param = true;
%!
%! sp = sp_evaluate_col (space, msh_col, 'value', true, 'gradient', true, 'hessian', true, 'third_derivative', true, 'fourth_derivative', true);
%!
%!
%! d2N11 = [0.573459956141965	-0.607197955911347	-0.567748083329291	0.471997865817610	0.179467816936778	0.0156200851868257	8.80693483837789	-4.22763190209158	-10.7500296797696	3.83603270738046	1.91421961937929	0.178161222783735	9.77086365868864	-2.11260801624193	-12.4441541341939	2.83578138437493	1.75909044710840	0.170522858772948	2.63086901080185	-0.0460475411773997	-3.39523352466546	0.511779941704709	0.406941558774687	0.0408936900280367	0.322957341734528	0.0517131332331411	-0.416401442808586	0.0383634194712236	0.0432685100506024	0.00451127157813988	0.0152815564594597	0.00498301046136600	-0.0194853228885160	0.000854274146194540	0.00177584536980883	0.000192578314334106];
%!
%! d2N12 = [0.580175813993755	1.22520762022496	-1.52774915235479	-0.886310589876210	-0.148729834526207	-0.00823562134564798	2.78377602674762	4.83186042158289	-4.30598256828398	-1.98738592059724	-0.227832782700092	-0.00598053864862797	0.199511177032473	-3.32094898240587	0.995910217716876	2.21885096402288	0.570624206068749	0.0427259781530717	-0.507396065893856	-2.84574117594319	0.974992102934703	1.32218376386310	0.294973967143568	0.0201429666577987	-0.121692447202472	-0.580155426319084	0.156705911308667	0.226767999157541	0.0494762036107502	0.00330073599062871	-0.00830210190433342	-0.0382900007542971	0.00750969549811197	0.0130068999081149	0.00284230653475599	0.000188230604878068];
%!
%!
%! d2N22 = [0.380134909284880	3.62536282876746	3.14018691982852	0.688995126876008	0.0682147453807321	0.00255069824111406	0.307301802306890	-0.336947857425302	-1.70100879154064	-0.580150760909668	-0.0700826041118709	-0.00270056162569342	-0.609501166415254	-5.46712330510247	-3.38750105526799	-0.265086537526570	0.0477403240138286	0.00600517689818393	-0.137985869147205	0.249124771636167	1.59878855458837	0.719388564251251	0.120186549884997	0.00701978080156724	0.00505302910668552	0.460092597202850	0.717761270774944	0.237466235435763	0.0336857618042125	0.00176715426135508	0.00193862496731436	0.0499100443413345	0.0663717636311214	0.0202027789781010	0.00270296387006736	0.000135531938941016];
%!
%! assert(d2N11(:)', permute(sp.shape_function_hessians(1,1,1,:,1),[1,4,3,2,5]),5e-14)
%! assert(d2N12(:)', permute(sp.shape_function_hessians(2,1,1,:,1),[1,4,3,2,5]),5e-14)
%! assert(d2N22(:)', permute(sp.shape_function_hessians(2,2,1,:,1),[1,4,3,2,5]),5e-14)
%!
%! d3N111 = [-7.96875045496924	16.0163816315465	-9.14750886015214	-1.36194084518250	2.04474864222025	0.326286939139424	-149.293349899929	234.170132192381	-64.9698004340314	-39.4393214858016	16.3577572146734	3.30679661331500	-179.263036591341	246.943105223888	-34.2783780719956	-48.2246503217827	12.1695790118280	2.95212507915367	-51.0422512972795	62.9608634126590	-1.71956961625788	-13.3390665565494	2.25405360351921	0.666608520142586	-6.57178757486280	7.25240152810421	0.590434131042898	-1.63040271313550	0.180796908550610	0.0692344590093065	-0.324572541132607	0.318386112449157	0.0621720451020285	-0.0751761867896050	0.00492526253786174	0.00277491993042141];
%!
%! d3N222 = [-5.34120438261789	-41.5865732312151	-31.7214413000696	-6.22560420432887	-0.551758476328870	-0.0184447405635832	6.19474147358912	59.3356978013014	45.5695030807670	8.28817471229815	0.622180273777625	0.0157558740358351	4.05552062136566	-2.25597835011926	-16.5830373786032	-4.53199278168536	-0.258918955042395	0.0136622851540697	-1.42247185118262	-15.8054329151897	-7.16861728209418	1.27434477113935	0.613183120676638	0.0526684578703524	-0.325428081566110	0.159402228077032	3.84854151395730	1.94395088523077	0.349562991203043	0.0216500141173212	-0.00596744526614005	0.393186031917795	0.738839333839710	0.269054099746001	0.0409730536669652	0.00227875214178018];
%!
%! d3N112 = [-12.2324451106732	10.5730405871177	11.6434361346994	-6.90063710338542	-2.63555658866126	-0.218752973046683	-60.1988879422529	36.4603182119669	44.6328675004091	-15.1461226896624	-4.77828362504010	-0.250506760496169	0.684504444789785	19.0183531097843	-34.1690353608250	4.53115713939314	6.90250240858655	0.877011650148332	14.2756522129099	6.59865246697346	-27.1941997783871	2.85207913763608	3.66608559321156	0.422052354232268	3.33925674037609	1.54657262438814	-5.33349388190071	0.311449669063323	0.592127780272384	0.0680801092127060	0.229890514257855	0.130102140427238	-0.338726347249810	0.00542848427816578	0.0322399960544277	0.00378715139197541];
%!
%! d3N122 = [-10.2158348934760	-19.4703747726189	20.4544618092389	11.1475590988758	1.72422088595848	0.0872562807216223	-6.88707383626682	5.55750294736357	-3.42917913562327	-6.96553710190580	-1.53952008739371	-0.0899015045771770	17.5190668327793	33.3187197151097	-28.6493478924464	-10.4658394376916	-0.271671577383747	0.0929720712051276	3.33250484236564	-6.85767797904385	-3.25281679571102	5.27879436921648	1.86922711572357	0.163844934867458	-0.347930211134240	-4.94750355411288	0.392019286658680	2.15589121544971	0.570921193227959	0.0433502732354826	-0.0744120832794841	-0.525095293482902	0.0446935633967821	0.186976973131795	0.0463770467883126	0.00335570083344077];
%!
%! assert(d3N111(:)',permute(sp.shape_function_third_derivatives(1,1,1,1,:,1),[1,5,3,4,2,6]),5e-12)
%! assert(d3N222(:)',permute(sp.shape_function_third_derivatives(2,2,2,1,:,1),[1,5,3,4,2,6]),5e-12)
%! assert(d3N112(:)',permute(sp.shape_function_third_derivatives(1,1,2,1,:,1),[1,5,3,4,2,6]),5e-12)
%! assert(d3N112(:)',permute(sp.shape_function_third_derivatives(1,2,1,1,:,1),[1,5,3,4,2,6]),5e-12)
%! assert(d3N112(:)',permute(sp.shape_function_third_derivatives(2,1,1,1,:,1),[1,5,3,4,2,6]),5e-12)
%! assert(d3N122(:)',permute(sp.shape_function_third_derivatives(1,2,2,1,:,1),[1,5,3,4,2,6]),5e-12)
%! assert(d3N122(:)',permute(sp.shape_function_third_derivatives(2,1,2,1,:,1),[1,5,3,4,2,6]),5e-12)
%! assert(d3N122(:)',permute(sp.shape_function_third_derivatives(2,2,1,1,:,1),[1,5,3,4,2,6]),5e-12)
%!
%! d4N1111 = [44.1180280950069	-108.244047615425	111.142556453888	-59.1209268658035	7.03753276970412	5.06778294876522	1575.45971900995	-3037.44324182025	1962.37614998410	-511.003083965015	-31.2728205499819	42.1856430106908	2245.16945132736	-3992.77263303888	2147.26362622167	-357.762538023189	-75.8198835234176	32.9591460509941	711.486833906808	-1189.46904346949	551.651736818827	-54.3507762117772	-25.5556693063165	6.54532675476186	99.6914766711146	-157.233688597538	62.9199136553354	-2.21346845021195	-3.44329531028735	0.588372759243432	5.28975527001960	-7.87312592174654	2.69007470104311	0.0827581310309550	-0.167460966547916	0.0198190955688997];
%!
%! d4N2222 = [51.4820385685158	322.027336608727	220.075703437930	40.8168481706751	3.61861797055599	0.130201166680407	-115.470999068252	-620.147353565049	-348.871303078342	-52.0937561650415	-3.87331958556005	-0.136080483145630	45.1838767377175	376.806698016322	188.591327246021	10.5540025491838	-1.03250316726313	-0.00425764086142088	17.4895152154742	-28.0001742547625	-81.6429760191845	-13.0471824013192	1.74780019946677	0.339609366379946	-2.06021542608660	-30.5489752894562	-4.01615050342946	9.90332651689110	2.93347791835102	0.234403345269443	-0.331593083506372	0.406309455664727	5.49577410343987	2.86920779283453	0.536279336790945	0.0344860083678705];
%!
%! d4N1112 = [179.628145797838	-325.215751047759	143.733229239584	35.2197270819888	-27.9023747598647	-4.49314780974755	977.735575913703	-1434.28146668240	359.343109252282	147.349471606826	-47.7607081649235	-6.13933017631656	-94.6917648306474	123.350473027718	82.4660313221173	-155.267259271555	35.0497796418894	13.2989285821947	-303.404680092580	351.903890604351	40.8567420106259	-112.264491526288	16.5184755900563	6.37363919328932	-71.5888119172864	73.9394703085844	14.8677949836918	-21.3974334830130	1.98721016686728	0.984206313294523	-5.07048288733178	4.62547522420046	1.45134955670501	-1.32252923304646	0.0657086562243000	0.0517978087229393];
%!
%! d4N1122 = [216.190315728420	-156.989972558060	-184.911035828387	82.5427541482144	30.8445101612578	2.38305101512650	126.390265554431	-105.692307961913	55.4697538438055	-30.1994815262620	-22.5426264738965	-2.23790114026812	-397.168639961636	213.612382766665	295.527116395474	-100.429973480644	-21.2639128472052	0.499883907006400	-59.7165753728081	118.130119056785	-70.5445207835155	-14.2450831510937	16.9132034768310	2.91383000332956	13.8850215381681	29.3739468871348	-45.0998573364712	-1.15385117664986	5.84930104528264	0.813911596500174	2.30005670287112	2.74647810292995	-4.55764418296842	-0.165704968084672	0.470175192864634	0.0630116267657426];
%!
%! d4N1222 = [141.346882985822	230.971013464312	-214.234354678741	-106.070110560986	-14.8317095289673	-0.676410002059919	-171.314203660618	-335.680245148067	315.745065011419	151.784744410136	19.1428404289179	0.721448273436038	-89.0250771123732	99.3648590893468	4.80708334896772	-58.1782618832832	-9.97786007000112	-0.132358035081743	46.0299326638353	106.957061045894	-91.3979897722300	-22.4880279150468	4.68770579558457	0.943406015188789	8.60885225347875	-14.4120604149283	-15.6887569138554	10.8797308314467	4.82784667294630	0.466259146469634	-0.0177714706950790	-5.00050168250458	-0.818565762407857	1.98588932548802	0.621950556426653	0.0516932927313724];
%!
%! assert(d4N1111(:)',permute(sp.shape_function_fourth_derivatives(1,1,1,1,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N2222(:)',permute(sp.shape_function_fourth_derivatives(2,2,2,2,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1112(:)',permute(sp.shape_function_fourth_derivatives(1,1,1,2,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1112(:)',permute(sp.shape_function_fourth_derivatives(1,1,2,1,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1112(:)',permute(sp.shape_function_fourth_derivatives(1,2,1,1,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1112(:)',permute(sp.shape_function_fourth_derivatives(2,1,1,1,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1122(:)',permute(sp.shape_function_fourth_derivatives(1,1,2,2,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1122(:)',permute(sp.shape_function_fourth_derivatives(1,2,1,2,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1122(:)',permute(sp.shape_function_fourth_derivatives(1,2,2,1,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1122(:)',permute(sp.shape_function_fourth_derivatives(2,2,1,1,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1122(:)',permute(sp.shape_function_fourth_derivatives(2,2,1,1,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1122(:)',permute(sp.shape_function_fourth_derivatives(2,1,1,2,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1222(:)',permute(sp.shape_function_fourth_derivatives(1,2,2,2,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1222(:)',permute(sp.shape_function_fourth_derivatives(2,1,2,2,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1222(:)',permute(sp.shape_function_fourth_derivatives(2,2,1,2,1,:,1),[1,6,3,4,5,2,7]),5e-11)
%! assert(d4N1222(:)',permute(sp.shape_function_fourth_derivatives(2,2,2,1,1,:,1),[1,6,3,4,5,2,7]),5e-11)
