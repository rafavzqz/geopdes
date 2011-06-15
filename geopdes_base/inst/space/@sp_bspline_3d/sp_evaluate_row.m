function [sp, elem_list] = sp_evaluate_row (space, msh, rownum, varargin)

value = true;
gradient = true;
if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_row: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'value'))
      value = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    else
      error ('sp_evaluate_row: unknown option %s', varargin {ii});
    end
  end
end

nel_row = msh.nelu * msh.nelv;
elem_list = nel_row * (rownum-1) + (1:nel_row);

spu = space.spu;
spv = space.spv;
spw = space.spw;

nsh  = spu.nsh' * spv.nsh * spw.nsh(rownum);
nsh  = nsh(:)';
ndof = spu.ndof * spv.ndof * spw.ndof;
ndof_dir = [spu.ndof, spv.ndof spw.ndof];

connectivity = space.connectivity(:,elem_list);

shp_u = reshape (spu.shape_functions, msh.nqnu, 1, 1, spu.nsh_max, 1, 1, msh.nelu, 1);
shp_u = repmat  (shp_u, [1, msh.nqnv, msh.nqnw, 1, spv.nsh_max, spw.nsh_max, 1, msh.nelv]);
shp_u = reshape (shp_u, msh.nqn, space.nsh_max, nel_row);

shp_v = reshape (spv.shape_functions, 1, msh.nqnv, 1, 1, spv.nsh_max, 1, 1, msh.nelv);
shp_v = repmat  (shp_v, [msh.nqnu, 1, msh.nqnw, spu.nsh_max, 1, spw.nsh_max, msh.nelu, 1]);
shp_v = reshape (shp_v, msh.nqn, space.nsh_max, nel_row);

shp_w = reshape (spw.shape_functions(:,:,rownum), 1, 1, msh.nqnw, 1, 1, spw.nsh_max, 1, 1);
shp_w = repmat  (shp_w, [msh.nqnu, msh.nqnv, 1, spu.nsh_max, spv.nsh_max, 1, msh.nelu, msh.nelv]);
shp_w = reshape (shp_w, msh.nqn, space.nsh_max, nel_row);

sp = struct('nsh_max', space.nsh_max, 'nsh', nsh, 'ndof', ndof,  ...
            'ndof_dir', ndof_dir, 'connectivity', connectivity, ...
            'ncomp', 1);
if (value)
  sp.shape_functions = shp_u .* shp_v .* shp_w;
end

if (gradient)
  shg_u = reshape (spu.shape_function_gradients, ...
            msh.nqnu, 1, 1, spu.nsh_max, 1, 1, msh.nelu, 1);
  shg_u = repmat  (shg_u, [1, msh.nqnv, msh.nqnw, 1, spv.nsh_max, spw.nsh_max, 1, msh.nelv]);
  shg_u = reshape (shg_u, msh.nqn, sp.nsh_max, nel_row);

  shg_v = reshape (spv.shape_function_gradients, ...
            1, msh.nqnv, 1, 1, spv.nsh_max, 1, 1, msh.nelv);
  shg_v = repmat  (shg_v, [msh.nqnu, 1, msh.nqnw, spu.nsh_max, 1, spw.nsh_max, msh.nelu, 1]);
  shg_v = reshape (shg_v, msh.nqn, sp.nsh_max, nel_row);

  shg_w = reshape (spw.shape_function_gradients(:,:,rownum), ...
            1, 1, msh.nqnw, 1, 1, spw.nsh_max, 1, 1);
  shg_w = repmat  (shg_w, [msh.nqnu, msh.nqnv, 1, spu.nsh_max, spv.nsh_max, 1, msh.nelu, msh.nelv]);
  shg_w = reshape (shg_w, msh.nqn, sp.nsh_max, nel_row);

  shape_fun_grads(1,:,:,:) = shg_u .* shp_v .* shp_w ;
  shape_fun_grads(2,:,:,:) = shp_u .* shg_v .* shp_w ;
  shape_fun_grads(3,:,:,:) = shp_u .* shp_v .* shg_w ;
  clear  shg_u shg_v shg_w

  JinvT = geopdes_invT__ (msh.geo_map_jac(:,:,:,elem_list));
  JinvT = reshape (JinvT, [3, 3, msh.nqn, nel_row]);
  shape_fun_grads = reshape (shape_fun_grads, [3, msh.nqn, sp.nsh_max, nel_row]);
  sp.shape_function_gradients = geopdes_prod__ (JinvT, shape_fun_grads);

  clear shape_fun_grads
end

clear shp_u shp_v shp_w

end
