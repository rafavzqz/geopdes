% SP_EVALUATE_COL_PARAM: compute the basis functions, in the parametric domain, in one column of the mesh.
%
%     sp = sp_evaluate_col_param (space, msh_col, 'option1', value1, ...)
%
% INPUTS:
%
%    space:   object defining the space of discrete functions (see sp_scalar)
%    msh_col: msh structure containing (in the field msh.qn) the points 
%              along each parametric direction in the parametric 
%              domain at which to evaluate, i.e. quadrature points 
%              or points for visualization (see msh_cartesian/msh_evaluate_col)
%   'option', value: additional optional parameters, currently available options are:
%            
%              Name             |   Default value |  Meaning
%           --------------------+-----------------+----------------------------------
%            value              |      true       |  compute shape_functions
%            gradient           |      false      |  compute shape_function_gradients
%            hessian            |      false      |  compute shape_function_hessians
%            third_derivative   |      false      |  compute shape_function_third_derivatives
%            fourth_derivative  |      false      |  compute shape_function_fourth_derivatives
%
% OUTPUT:
%
%    sp: struct representing the discrete function space, with the following fields:
%              (see the article for a detailed description)
%
%    FIELD_NAME      (SIZE)                                 DESCRIPTION
%    ncomp           (scalar)                               number of components of the functions of the space (actually, 1)
%    ndof            (scalar)                               total number of degrees of freedom
%    ndof_dir        (1 x ndim vector)                      degrees of freedom along each direction
%    nsh_max         (scalar)                               maximum number of shape functions per element
%    nsh             (1 x msh_col.nel vector)               actual number of shape functions per each element
%    connectivity    (nsh_max x msh_col.nel vector)         indices of basis functions that do not vanish in each element
%    shape_functions (msh_col.nqn x nsh_max x msh_col.nel)  basis functions evaluated at each quadrature node in each element
%    shape_function_gradients
%         (ndim x msh_col.nqn x nsh_max x msh_col.nel)      basis function gradients evaluated at each quadrature node in each element
%    shape_function_hessians
%         (ndim x ndim x msh_col.nqn x nsh_max x msh_col.nel)  basis function hessians evaluated at each quadrature node in each element
%    shape_function_third_derivatives
%       (ndim x ndim x ndim x msh_col.nqn x nsh_max x msh_col.nel) basis function third derivatives evaluated at each quadrature node in each element
%    shape_function_fourth_derivatives
%       (ndim x ndim x ndim x rdim x msh_col.nqn x nsh_max x msh_col.nel) basis function fourth derivatives evaluated at each quadrature node in each element
%
% Copyright (C) 2009, 2010, 2011 Carlo de Falco
% Copyright (C) 2011, 2015, 2019 Rafael Vazquez
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

function sp = sp_evaluate_col_param (space, msh, varargin)

value = true;
gradient = false;
hessian = false;
third_derivative = false;
fourth_derivative = false;

if (~isempty (varargin))
  if (~rem (length (varargin), 2) == 0)
    error ('sp_evaluate_col_param: options must be passed in the [option, value] format');
  end
  for ii=1:2:length(varargin)-1
    if (strcmpi (varargin {ii}, 'value'))
      value = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'gradient'))
      gradient = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'hessian'))
      hessian = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'third_derivative'))
      third_derivative = varargin {ii+1};
    elseif (strcmpi (varargin {ii}, 'fourth_derivative'))
      fourth_derivative = varargin {ii+1};      
    else
      error ('sp_evaluate_col_param: unknown option %s', varargin {ii});
    end
  end
end

sp_univ = space.sp_univ;

elem_list{1} = msh.colnum;
for idim = 2:msh.ndim
  elem_list{idim} = 1:msh.nel_dir(idim);
end

for idim = 1:msh.ndim
  nsh_dim{idim} = sp_univ(idim).nsh(elem_list{idim});
end

[nsh_grid{1:msh.ndim}] = ndgrid (nsh_dim{:});
nsh = 1;
for idim = 1:msh.ndim
  nsh = nsh .* nsh_grid{idim};
end
nsh = nsh(:)';

for idim = 1:msh.ndim
  csize = ones (1, 2*msh.ndim);
  csize([idim, msh.ndim+idim]) = [sp_univ(idim).nsh_max, msh.nel_dir(idim)];
  crep = [sp_univ.nsh_max, msh.nel_dir];
  crep([idim, msh.ndim+idim]) = 1;

  conn{idim} = reshape (sp_univ(idim).connectivity(:,elem_list{idim}), csize);
  conn{idim} = repmat (conn{idim}, crep);
  conn{idim} = reshape (conn{idim}, [], msh.nel);
end

connectivity = zeros (space.nsh_max, msh.nel);
indices = ones (size (conn{1}));
for idim = 1:msh.ndim
  indices = indices & conn{idim} ~= 0;
end
for idim = 1:msh.ndim
  conn{idim} = conn{idim}(indices);
end
connectivity(indices) = sub2ind ([space.ndof_dir, 1], conn{:}); % The extra 1 makes things work in any dimension
connectivity = reshape (connectivity, space.nsh_max, msh.nel);

clear conn csize crep indices

sp = struct('nsh_max', space.nsh_max, 'nsh', nsh, 'ndof', space.ndof,  ...
            'ndof_dir', space.ndof_dir, 'connectivity', connectivity, ...
            'ncomp', 1, 'degree', space.degree);

shp = cell(1,msh.ndim); shg = cell(1,msh.ndim); shh = cell(1,msh.ndim); shtd = cell(1,msh.ndim); shfd = cell(1,msh.ndim);
for idim = 1:msh.ndim
  ssize = ones (1, 3*msh.ndim);
  ssize([idim, msh.ndim+idim, 2*msh.ndim+idim]) = [msh.nqn_dir(idim), sp_univ(idim).nsh_max, msh.nel_dir(idim)];
  srep = [msh.nqn_dir, sp_univ.nsh_max, msh.nel_dir];
  srep([idim, msh.ndim+idim, 2*msh.ndim+idim]) = 1;
  shp{idim} = reshape (sp_univ(idim).shape_functions(:,:,elem_list{idim}), ssize);
  shp{idim} = repmat (shp{idim}, srep);
  shp{idim} = reshape (shp{idim}, msh.nqn, space.nsh_max, msh.nel);  
  shg{idim} = reshape (sp_univ(idim).shape_function_gradients(:,:,elem_list{idim}), ssize);
  shg{idim} = repmat (shg{idim}, srep);
  shg{idim} = reshape (shg{idim}, msh.nqn, space.nsh_max, msh.nel);  
  shh{idim} = reshape (sp_univ(idim).shape_function_hessians(:,:,elem_list{idim}), ssize);
  shh{idim} = repmat (shh{idim}, srep);
  shh{idim} = reshape (shh{idim}, msh.nqn, space.nsh_max, msh.nel);
  
  shtd{idim} = reshape (sp_univ(idim).shape_function_third_derivatives(:,:,elem_list{idim}), ssize);
  shtd{idim} = repmat (shtd{idim}, srep);
  shtd{idim} = reshape (shtd{idim}, msh.nqn, space.nsh_max, msh.nel);  
  shfd{idim} = reshape (sp_univ(idim).shape_function_fourth_derivatives(:,:,elem_list{idim}), ssize);
  shfd{idim} = repmat (shfd{idim}, srep);
  shfd{idim} = reshape (shfd{idim}, msh.nqn, space.nsh_max, msh.nel);  

end

if (value)
  sp.shape_functions = 1;
  for idim = 1:msh.ndim
    sp.shape_functions = sp.shape_functions .* shp{idim};
  end
end

if (gradient)
  for idim = 1:msh.ndim
    shape_fun_grad = shg{idim};
    for jdim = setdiff (1:msh.ndim, idim)
      shape_fun_grad = shape_fun_grad .* shp{jdim};
    end
    sp.shape_function_gradients(idim,:,:,:) = shape_fun_grad;
  end
  sp.shape_function_gradients = reshape (sp.shape_function_gradients, ...
                                msh.ndim, msh.nqn, sp.nsh_max, msh.nel);
end

if (hessian && isfield (msh, 'geo_map_der2'))
  for idim = 1:msh.ndim
    shape_fun_hess = shh{idim};
    for jdim = setdiff (1:msh.ndim, idim)
      shape_fun_hess = shape_fun_hess .* shp{jdim};
    end
    sp.shape_function_hessians(idim,idim,:,:,:) = shape_fun_hess;
    
    for jdim = setdiff (1:msh.ndim, idim)
      shape_fun_hess = shg{idim} .* shg{jdim};
      for kdim = setdiff (1:msh.ndim, [idim, jdim])
        shape_fun_hess = shape_fun_hess .* shp{kdim};
      end
      sp.shape_function_hessians(idim,jdim,:,:,:) = shape_fun_hess;
    end
  end
end

if (third_derivative && isfield (msh, 'geo_map_der3'))
  for idim = 1:msh.ndim
    shape_fun_third = shtd{idim};
    for jdim = setdiff (1:msh.ndim, idim)
      shape_fun_third = shape_fun_third .* shp{jdim};
    end
    sp.shape_function_third_derivatives(idim,idim,idim,:,:,:) = shape_fun_third;
    
    for jdim = setdiff (1:msh.ndim, idim)
      shape_fun_third = shh{idim} .* shg{jdim};
      for kdim = setdiff (1:msh.ndim, [idim, jdim])
        shape_fun_third = shape_fun_third .* shp{kdim};
      end
      sp.shape_function_third_derivatives(idim,idim,jdim,:,:,:) = shape_fun_third;
      sp.shape_function_third_derivatives(idim,jdim,idim,:,:,:) = shape_fun_third;
      sp.shape_function_third_derivatives(jdim,idim,idim,:,:,:) = shape_fun_third;
    end
  end
end

if (fourth_derivative && isfield (msh, 'geo_map_der4'))    
  for idim = 1:msh.ndim
    shape_fun_fourth = shfd{idim};
    for jdim = setdiff (1:msh.ndim, idim)
      shape_fun_fourth = shape_fun_fourth .* shp{jdim};
    end
    sp.shape_function_fourth_derivatives(idim,idim,idim,idim,:,:,:) = shape_fun_fourth;
    
    for jdim = setdiff (1:msh.ndim, idim)
      shape_fun_fourth = shh{idim} .* shh{jdim};
      for kdim = setdiff (1:msh.ndim, [idim, jdim])
        shape_fun_fourth = shape_fun_fourth .* shp{kdim};
      end
      sp.shape_function_fourth_derivatives(idim,idim,jdim,jdim,:,:,:) = shape_fun_fourth;
      sp.shape_function_fourth_derivatives(idim,jdim,idim,jdim,:,:,:) = shape_fun_fourth;
      sp.shape_function_fourth_derivatives(jdim,idim,idim,jdim,:,:,:) = shape_fun_fourth;
      sp.shape_function_fourth_derivatives(jdim,idim,jdim,idim,:,:,:) = shape_fun_fourth;
    end
    
    for jdim = setdiff (1:msh.ndim, idim)
      shape_fun_fourth = shtd{idim} .* shg{jdim};
      for kdim = setdiff (1:msh.ndim, [idim, jdim])
        shape_fun_fourth = shape_fun_fourth .* shp{kdim};
      end
      sp.shape_function_fourth_derivatives(idim,idim,idim,jdim,:,:,:) = shape_fun_fourth;
      sp.shape_function_fourth_derivatives(idim,idim,jdim,idim,:,:,:) = shape_fun_fourth;
      sp.shape_function_fourth_derivatives(idim,jdim,idim,idim,:,:,:) = shape_fun_fourth;
      sp.shape_function_fourth_derivatives(jdim,idim,idim,idim,:,:,:) = shape_fun_fourth;
      
    end
  end
end

clear shp shg shh shtd shfd

if (value || gradient || hessian)
  if (strcmpi (space.space_type, 'NURBS'))
    sp = bsp_2_nrb__ (sp, msh, space.weights);
  end
end

%test = scriptFourthDerivative_Ale(space, 1, space.nsh_max, sp);

end


%!test
%! base = 1; height = 1;
%! p11 =[0 0]; p12 =[base 0]; p21 =[0 height]; p22 =[base height];
%! srf = nrb4surf(p11,p12,p21,p22);
%! degree     = [5 5];       % Degree of the splines
%! regularity = [4 4];       % Regularity of the splines
%! nsub       = [3 3];       % Number of subdivisions
%! nquad      = [1 1];       % Points for the Gaussian quadrature rule
%! geo_name = 'geo_plate_with_hole.txt';
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
%! sp = sp_evaluate_col_param (space, msh_col, 'value', value_param, 'gradient', grad_param, 'hessian', hessian_param, ...
%!                            'third_derivative', third_param, 'fourth_derivative', fourth_param);
%!
%!
%! assert (msh_col.nel, 3)
%! assert (sp.ndof, 120)
%! % check values of shapes and their derivatives at one Gauss point in a element
%!  N = [0.00101438869526536	0.0131160438467795	0.0132316455223460	0.00326140938786764	0.000368329784979418	1.60951130638220e-05 ...
%!        	0.0133772509188119	0.173197472780233	0.175173593410800	0.0433258041472770	0.00490511607342200	0.000214619141218335	... 
%!          0.0140065476093932	0.181825975886514	0.184840795672781	0.0460256131519473	0.00523589224041652	0.000229667490082379	...
%!          0.00362080409282217	0.0471899354948230	0.0483360115804357	0.0121545553937682	0.00139230697163462	6.12911382029609e-05	...
%!          0.000422661956360565	0.00553032003749955	0.00570691115607263	0.00144876925451682	0.000167053786586711	7.37876280353217e-06	...
%!          1.87849758382473e-05	0.000246759435266767	0.000256510503555254	6.57208361450235e-05	7.62584300156177e-06	3.37907468210412e-07];
%!
%! assert(N(:)',sp.shape_functions(1,:,1),1e-15)
%!
%! dN1 = [-0.0604346887411712	-0.177212131495881	0.127104825208900	0.0904640867772802	0.0163621493905147	0.000972507821955720	-0.796982457774195 ...
%!	-2.34008773374306	1.68273922801908	1.20175937453755	0.217897778686647	0.0129678364326250	-0.834474348011773	-2.45666827014150	1.77560356992523 ...
%!	1.27664594259280	0.232591700492315	0.0138770960892531	-0.215718263978903	-0.637587762886601	0.464321712129524	0.337139752516628	0.0618498302237861 ...
%!	0.00370336704579782	-0.0251811203088213	-0.0747206865994362	0.0548212951858705	0.0401855676409133	0.00742095568658847	0.000445843686486470 ...
%!	-0.00111916090261428	-0.00333399049295498	0.00246407165787319	0.00182294668256474	0.000338759414819350	2.04172319031752e-05];
%!
%!
%! dN2 = [-0.0304788104851832	-0.394090959991200	-0.397564383519881	-0.0979938746472445	-0.0110670138230388	-0.000483601506108565	-0.0938187470865016 ...
%!	-1.21468678380992	-1.22854592144416	-0.303857099360868	-0.0344010774048622	-0.00150518959777974	0.0636640250854994	0.826454442083319	...
%!  0.840157716284336	0.209200430541214	0.0237987250131006	0.00104390869596750	0.0492831509192391	0.642307248125047	0.657906612035443	0.165436950597467	...
%!  0.0189508387777741	0.000834240223064400	0.00927891736192305	0.121409986966438	0.125286783473506	0.0318055836037640	0.00366741850646259	...
%! 0.000161989810667509	0.000562676133902558	0.00739131347495188	0.00768339228587134	0.00196857032542106	0.000228420834544284	1.01215178271518e-05];
%!
%! assert(dN1(:)', permute(sp.shape_function_gradients(1,1,:,1),[1,3,2,4]),5e-15)
%! assert(dN2(:)', permute(sp.shape_function_gradients(2,1,:,1),[1,3,2,4]),5e-15)
%!
%! d2N11 = [2.86933700654376	-1.05818773954205	-3.34937922321075	0.986846497528920	0.508489448017523	0.0471594753054431	37.8393817737958 ...
%!	-13.9733782805790	-44.3423906145631	13.1096446322916	6.77164830635569	0.628844671685583	39.6194334352628	-14.6695162551337	...
%!  -46.7894880937217	13.9265604939332	7.22829394684277	0.672937076245122	10.2419390372465	-3.80723118599998	-12.2354874642743	...
%!  3.67775982493375	1.92211825475096	0.179586058641744	1.19555708605990	-0.446180031706946	-1.44461319059507	0.438372707782715	...
%!  0.230622369395026	0.0216201390346285	0.0531358704915511	-0.0199082750916335	-0.0649315272006293	0.0198859968962046	0.0105276870257699	0.000990085552018459];
%!
%! d2N12 = [1.81531543386700	5.31769263623823	-3.82602185714106	-2.71984524706367	-0.491819008823616	-0.0292289149918344	5.58243557899033	...
%! 16.3205262748091	-11.8938289052578	-8.45112548200957	-1.53076722911548	-0.0910604390326279	-3.80031738851064	-11.2620740600474	7.97330262631940	...
%! 5.77850277406712	1.05444235652424	0.0629546616658333	-2.93807140516835	-8.70313022114258	6.29447359587945	4.58244329565577	0.841111305854630	...
%! 0.0503746417834245	-0.553036825280792	-1.64329486853791	1.20051450078485	0.881451556762397	0.162828105003659	0.00978395149331023	...
%! -0.0335327016778383	-0.0999947173783087	0.0736725124840521	0.0545690529640001	0.0101430200627028	0.000611389969953273];
%!
%!
%! d2N22 = [0.733193221331904	9.48018691852152	9.56374302013780	2.35732443227520	0.266226253133570	0.0116334378035215	-0.219515910202459	...
%! -2.84210867485882	-2.87453610868073	-0.710960973249628	-0.0804912030162901	-0.00352182346113885	-0.917177491207438	-11.9063381677366	...
%! -12.1037547564009	-3.01385163419410	-0.342856972554550	-0.0150390987301111	0.249015730420845	3.24542172249833	3.32424150016663	...
%! 0.835912524325701	0.0957539619994726	0.00421521218952204	0.141151833799527	1.84690105897632	1.90587528138785	0.483829769748594	...
%! 0.0557891429900861	0.00246420546069162	0.0134728752975037	0.176979684640456	0.183973301677788	0.0471360004640916	0.00546937258174133	0.000242352464038030];
%!
%! assert(d2N11(:)', permute(sp.shape_function_hessians(1,1,1,:,1),[1,4,3,2,5]),5e-14)
%! assert(d2N12(:)', permute(sp.shape_function_hessians(2,1,1,:,1),[1,4,3,2,5]),5e-14)
%! assert(d2N22(:)', permute(sp.shape_function_hessians(2,2,1,:,1),[1,4,3,2,5]),5e-14)
%!
%! d3N111 = [-101.350613763121	147.430156064144	-32.0007168048603	-24.0231347679074	8.23249772817786	1.72556471432208	-1336.56121900116	...
%! 1946.81648980418	-423.657098806170	-319.132469486911	109.633699411933	23.0094200417207	-1399.43613679865	2043.80469557870	-447.037214405434	...
%! -339.018925878456	117.026840434916	24.6227606691295	-361.765385237807	530.435826224148	-116.900578650608	-89.5289390392650	31.1192970222288	...
%! 6.57105203078037	-42.2294224013004	62.1632525583099	-13.8021569144635	-10.6714536293199	3.73380045449068	0.791080663965454	-1.87686321783557	...
%! 2.77368560799772	-0.620370306012279	-0.484091481935579	0.170444361944029	0.0362272201218906];
%!
%! d3N222 = [-13.2487157763057	-171.305869088070	-172.815718073961	-42.5965768464429	-4.81067726396889	-0.210214860797815	18.9286712370784	...
%! 245.072626748713	247.868816935160	61.3055632852474	6.94068834448987	0.303683857764312	-3.49347413241264	-45.3505289864747	-46.1024769490221	...
%! -11.4795803690361	-1.30592167406966	-0.0572830263412563	-3.50393271036684	-45.6667509043639	-46.7758342411035	-11.7622337835447	-1.34736644560623	...
%! -0.0593127986213632	1.07585722824619	14.0770530617214	14.5265537288665	3.68774348170359	0.425224038029787	0.0187821382507444	0.241574595635586	...
%! 3.17332379381950	3.29872243149090	0.845169275344031	0.0980682623893451	0.00434548655787843];
%!
%! d3N112 = [-86.1489861965118	31.9922360131239	100.514001638292	-29.7438890038355	-15.2952528263121	-1.41798753716397	-264.528639728833	...
%! 100.606481845794	309.357580939923	-93.1724033713140	-47.7171563412657	-4.42376443912575	180.972905713070	-63.9406346466879	-214.391531610059	...
%! 61.9933232847200	32.6140716105913	3.04427381989899	139.634330012742	-51.1103425713403	-166.988060332177	49.7131898271135	26.0981438195995	...
%! 2.44051404718372	26.2735575599015	-9.71198131924695	-31.7674188030308	9.58267832741669	5.05529212269552	0.474174524277660	1.59280063326137	...
%! -0.592608757657191	-1.94731334439780	0.593789105265513	0.314990668473306	0.0296353156158172];
%!
%! d3N122 = [-43.6496239902832	-127.671937471365	92.2899451444605	65.4902113724882	11.8381259567232	0.703432230080403	13.1776819586829	...
%! 39.6880455815338	-26.3104121121722	-19.3982351610114	-3.53914617659307	-0.211201476872619	54.5767149704975	160.005990319521	-117.146145637724	...
%! -83.8155170452038	-15.2553954410731	-0.909789156831090	-14.8874594827566	-44.5235067035239	31.2423977925847	23.0126367571585	4.23374064999226	...
%! 0.253818087237993	-8.41921969331642	-25.0812815826582	18.1763752146629	13.3869052607431	2.47444105071116	0.148723304884823	-0.803271322119616	...
%! -2.39896355195272	1.75918912483236	1.30537507332142	0.242723246287450	0.0146329090534230];
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
%! d4N1111 = [2335.36042349566	-4230.01520559041	2422.52568737248	-520.272158366309	-50.2808630942555	42.6794076148571	30797.5655848491	...
%! -55857.3874857932	32071.7879775671	-6911.49345449012	-669.599581155092	569.105527478213	32246.3539957214	-58640.1396454711	33841.7149125898	...
%! -7342.17703051152	-714.753983122404	609.009230704923	8335.93928942202	-15219.0818379237	8849.63472464443	-1938.93989274163	-190.064445181509	...
%! 162.525697096379	973.066843123193	-1783.56660923701	1044.85408468585	-231.113061066838	-22.8045868547197	19.5662636314847	43.2474152499197	...
%! -79.5816311301536	46.9634168247637	-10.4840322708354	-1.04100722661145	0.896029155848123];
%!
%! d4N2222 = [160.220975043871	2071.65689410552	2089.91598281819	515.134084774767	58.1770652241587	2.54219582746073	-289.455166725741	...
%! -3747.62375800094	-3790.38279198054	-937.477957099627	-106.136245739762	-4.64390028122505	163.154587272964	2117.99101952672	2153.10900089143	...
%! 536.127111919795	60.9900356114867	2.67527056626566	-32.2147309416604	-419.854550577225	-430.051327809630	-108.140546046654	-12.3875231383887	...
%! -0.545314652684657	-4.58218610585446	-59.9556105190150	-61.8700798904968	-15.7064771237014	-1.81107272208572	-0.0799950501527926	2.87652509273153	...
%! 37.7860324934729	39.2792041032902	10.0637677638927	1.16773792716628	0.0517434421900925];
%!
%! d4N1112 = [3040.54959401812	-4428.48606147947	967.148544704528	720.482203678773	-248.119882840565	-51.9192325258939	9311.99884056601	...
%! -13636.8312650288	3045.89951968251	2220.52579859826	-779.046283733301	-162.333309402340	-6425.47989713952	9307.33866775112	-1953.13705503769	...
%! -1559.69167548047	521.086398474120	110.889331847164	-4940.73112166055	7224.39260867014	-1570.54437168850	-1223.53880175721	420.686545406370	...
%! 89.1647849359203	-929.034300693402	1365.23819306206	-300.573504346120	-234.866082213958	81.6243242744001	17.3339571494983	-56.3053191068358	...
%! 83.1055626389189	-18.4729433899164	-14.5270132443833	5.08962596968530	1.08361934132283];
%!
%! d4N1122 = [2070.05673779532	-776.718124652077	-2413.48744776753	718.864227163247	368.552320088938	34.1474285512547	-632.924644050500	...
%! 192.681863450182	750.840078553517	-197.702927018457	-107.930666558410	-10.1284806739208	-2586.33376855720	985.416431357321	3048.60937859305	...
%! -923.649846439828	-475.483684183656	-44.1948838397945	710.624313831397	-242.515061489624	-853.620416821718	243.581477744445	130.455914598652	...
%! 12.2463739510037	400.445093025756	-145.352665457803	-484.760932213381	144.597446796022	76.6821258008953	7.19993900190462	38.1813262550265	...
%! -14.0560911070556	-46.7123036716889	14.1509755454121	7.52965725813709	0.708835141158010];
%!
%! d4N1222 = [788.161650331741	2299.49235378646	-1675.25946494326	-1185.27152056410	-214.124861273261	-12.7201718532330	-1127.38892012784	...
%! -3306.88861662476	2385.41897994042	1701.55551360167	308.445194263800	18.3546963027006	209.590769689004	631.668875723252	-423.618657873354	...
%! -313.624838695205	-57.4671620999680	-3.43727145352199	208.369106965905	611.973763697739	-454.490821492699	-327.554260251829	-60.0019851861459	...
%! -3.59036910450154	-64.3185027051784	-193.096935034881	136.550740628699	101.529782458290	18.8019217998291	1.13099476588722	-14.4136038076010	...
%! -43.1536471759577	31.3983644490292	23.3688840721991	4.34783385861134	0.262183932064408];
%!
%! assert(d4N1111(:)',permute(sp.shape_function_fourth_derivatives(1,1,1,1,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)
%! assert(d4N2222(:)',permute(sp.shape_function_fourth_derivatives(2,2,2,2,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)
%! assert(d4N1112(:)',permute(sp.shape_function_fourth_derivatives(1,1,1,2,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)
%! assert(d4N1112(:)',permute(sp.shape_function_fourth_derivatives(1,1,2,1,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)
%! assert(d4N1112(:)',permute(sp.shape_function_fourth_derivatives(1,2,1,1,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)
%! assert(d4N1112(:)',permute(sp.shape_function_fourth_derivatives(2,1,1,1,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)
%! assert(d4N1122(:)',permute(sp.shape_function_fourth_derivatives(1,1,2,2,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)
%! assert(d4N1122(:)',permute(sp.shape_function_fourth_derivatives(1,2,1,2,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)
%! assert(d4N1122(:)',permute(sp.shape_function_fourth_derivatives(1,2,2,1,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)
%! assert(d4N1122(:)',permute(sp.shape_function_fourth_derivatives(2,2,1,1,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)
%! assert(d4N1222(:)',permute(sp.shape_function_fourth_derivatives(1,2,2,2,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)
%! assert(d4N1222(:)',permute(sp.shape_function_fourth_derivatives(2,1,2,2,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)
%! assert(d4N1222(:)',permute(sp.shape_function_fourth_derivatives(2,2,1,2,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)
%! assert(d4N1222(:)',permute(sp.shape_function_fourth_derivatives(2,2,2,1,1,:,1),[1,6,3,4,5,2,7]),7.5e-11)