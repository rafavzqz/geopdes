function diff_ops = op_geom_exterior (knots, degrees, grad_curl, coo)

% OP_GEOM_EXTERIOR: computes the exterior derivatives d_0,...,d_(n-1)
% starting from the H^1-conforming space of 0-forms which is deduced by the
% input knot vector and polynomial degrees 
%
%   diff_ops = op_geom_exterior (knots, degrees, [grad_curl]);
%
% INPUT:
%
%   knots:      initial knot vector
%   degrees:    polynomial degree(s) of the (multivariate) space of 0-forms
%   grad_curl:  can be either 'grad' or 'curl', a flag which denotes which 
%               sequence of operators is sought between the one starting with
%               the gradient and the one starting with the rotated gradient.
%               Only valid in the 2D case. Default value is 'grad'.
%   coo:        true if COO representation of sparse matrices is desired,
%               false if Matlab sparse matrix representation is desired.
%               Default value is false.
%               
%
% OUTPUT:
%
%   diff_ops:   cell-array containing the (sparse) matrix representation of the
%               differential operators. The k-th entry corresponds to the
%               exterior derivative d_(k-1).
% 
% Copyright (C) 2020-2023 Bernard Kapidani, Rafael Vazquez
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

  
  if nargin < 3
    grad_curl = 'grad';
    coo = false;
  else
    assert(ismember(grad_curl,{'grad','curl'}), ...
           '    op_geom_exterior: Unrecognized input parameter!');
    if nargin < 4
      coo = false;
    end
  end
  ndim = numel(knots);
  
  Id       = cell(ndim,1);
  Id_deriv = cell(ndim,1);
  D_univ   = cell(ndim,1);
  
  if strcmpi(grad_curl,'curl')
    assert(ndim == 2, strcat('    op_geom_exterior: div-curl short exact sequence',...
                             ' only available in two dimensions!'));
  end

  if ~iscell(knots)
    knots = {knots};
  end
  
  if numel(degrees) == 1
    degrees = repmat(degrees,1,ndim);
  end
  
  %% Build matrices for univariate differential operators
  for kdim = 1:ndim
    degree = degrees(kdim);
    dim_univ = numel(knots{kdim})-degree-1;
    dim_d_univ = dim_univ-1;
    d_knots = knots{kdim}(2:end-1);
    Id{kdim} = speye(dim_univ);
    Id_deriv{kdim} = speye(dim_univ-1);

    scaling = degree./(abs(d_knots((1:dim_d_univ)+degree) - d_knots(1:dim_d_univ)).');

    d_univ_rows = repmat((1:dim_d_univ).',1,2);
    d_univ_cols = [ d_univ_rows(:,1), d_univ_rows(:,2)+1];
    d_univ_vals = bsxfun (@times, scaling(:), [ -ones(dim_d_univ,1), ones(dim_d_univ,1)]);
    
    D_univ{kdim}   = sparse(d_univ_rows(:), ...
                         d_univ_cols(:), ...
                         d_univ_vals(:), ...
                         dim_d_univ,dim_univ,numel(d_univ_vals));
  end

  %% Assemble multivariate topologic differential operators
  diff_ops = cell (ndim, 1);
  switch ndim
    case 1
      diff_ops{1} = D_univ{1};

    case 2
      if strcmpi(grad_curl, 'grad')  
        diff_ops{1} = vertcat(kron(Id{2},D_univ{1}),kron(D_univ{2},Id{1}));
        diff_ops{2} = horzcat(kron(kron(-1,D_univ{2}),Id_deriv{1}),...
                              kron(kron( 1,Id_deriv{2}),D_univ{1}));
      else        
        diff_ops{1} = vertcat(kron(D_univ{2},Id{1}),-kron(Id{2},D_univ{1}));
        diff_ops{2} = horzcat(kron(Id_deriv{2},D_univ{1}),...
                              kron(D_univ{2},Id_deriv{1}));
        
      end
      
    case 3
      diff_ops{1} = vertcat(kron(kron(Id{3},Id{2}),D_univ{1}),...
                            kron(kron(Id{3},D_univ{2}),Id{1}),...
                            kron(kron(D_univ{3},Id{2}),Id{1}));
      
      univ_c = cell(ndim,1); C = cell(ndim);
      dim_imm = [ size(Id_deriv{3},2)*size(Id_deriv{2},2)*size(D_univ{1},2)
                  size(Id_deriv{3},2)*size(D_univ{2},2)*size(Id_deriv{1},2)
                  size(D_univ{3},2)*size(Id_deriv{2},2)*size(Id_deriv{1},2) ];
      dim_dom = [ size(Id{3},1)*size(Id{2},1)*size(D_univ{1},1)
                  size(Id{3},1)*size(D_univ{2},1)*size(Id{1},1)
                  size(D_univ{3},1)*size(Id{2},1)*size(Id{1},1) ];
      
      for kdim = 1:ndim        
        transv_dir = setdiff(1:ndim,kdim);
        C(kdim,:) = circshift({0, -1, 1},kdim-1,2);
        C{kdim,kdim} = spalloc(dim_imm(kdim),dim_dom(kdim),0);
        univ_c{kdim} = Id{kdim};
        
        for mm = transv_dir
          nn = setdiff(transv_dir,mm);
          univ_c{mm} = Id_deriv{mm};
          univ_c{nn} = D_univ{nn};
          
          for ll = ndim:-1:1
            C{kdim,mm} = kron(C{kdim,mm},univ_c{ll});
          end
        end
      end
      
      diff_ops{2} = vertcat(horzcat(C{1,:}),horzcat(C{2,:}),horzcat(C{3,:}));
      diff_ops{3} = horzcat(kron(kron(Id_deriv{3},Id_deriv{2}),D_univ{1}),...
                            kron(kron(Id_deriv{3},D_univ{2}),Id_deriv{1}),...
                            kron(kron(D_univ{3},Id_deriv{2}),Id_deriv{1}));

      clear C
      
    otherwise
      error('    op_geom_exterior: Not implemented for dimension greater than three');
  end
  
  if coo
    for idim = 1:ndim
      [rows, cols, vals] = find (diff_ops{idim});
      diff_ops{idim} = [rows, cols, vals];
    end
  end

end
