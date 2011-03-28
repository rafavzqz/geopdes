% Copyright (C) 2011 Rafael Vazquez
%
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 3 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
%
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.

function [knotsp, knotsv1, degreev1, knotsv2, degreev2, fun_space, der2] = ...
          space_set (element_name, geometry, nsub_p, degree_p, regularity_p)

  knotsp = kntrefine (geometry.nurbs.knots, nsub_p, degree_p, regularity_p);

  switch lower (element_name)
   case 'th'
    fun_space = 'sp_bspline_th_2d_rafa';

    degreev1 = degree_p + 1;
    degreev2 = degreev1;
    regularity_v = regularity_p;
    nsub_v = nsub_p;
    knotsv1 = kntrefine (geometry.nurbs.knots, nsub_v, degreev1, regularity_v);
    knotsv2 = knotsv1;

    der2 = false;
   case 'sg'
    fun_space = 'sp_bspline_th_2d_rafa';

    degreev1 = degree_p + 1;
    degreev2 = degreev1;
    regularity_v = regularity_p+1;
    nsub_v = 2*nsub_p;
    knotsv1 = kntrefine (geometry.nurbs.knots, nsub_v, degreev1, regularity_v);
    knotsv2 = knotsv1;

    der2 = false;
   case 'ndl'
    fun_space = 'sp_bspline_ndl_2d_rafa';

    knotsv1{1}  = [knotsp{1}(1) knotsp{1} knotsp{1}(end)];
    knotsv1{2}  = sort ([knotsp{2}, unique(knotsp{2})]);
    knotsv2{1}  = sort ([knotsp{1}, unique(knotsp{1})]);
    knotsv2{2}  = [knotsp{2}(1) knotsp{2} knotsp{2}(end)];
    degreev1 = degree_p+1;
    degreev2 = degree_p+1;

    der2 = true;
   case 'rt'
    fun_space = 'sp_bspline_rt_2d_rafa';

    knotsv1{1}  = [knotsp{1}(1) knotsp{1} knotsp{1}(end)];
    knotsv1{2}  = knotsp{2};
    knotsv2{1}  = knotsp{1};
    knotsv2{2}  = [knotsp{2}(1) knotsp{2} knotsp{2}(end)];
    degreev1 = [degree_p(1)+1 degree_p(2)];
    degreev2 = [degree_p(1) degree_p(2)+1];

    der2 = true;
   otherwise
    error('space_set: Unknown element type')
  end

end
