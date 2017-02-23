% COLLOCATION_CSP: create the collocation Clustered Superconvergent Points from a given knot vector
%
% USAGE:
%
%  col_pts = collocation_csp (knots, degree)
%
% INPUT:
%
%   knots:  knots along each direction (in parametric space)
%   degree: degree of the B-splines used (from 3 to 7)
%
% OUTPUT:
%
%   col_pts: CSP collocation points for degree p in each parametric
%                   direction and continuity C^(p-1)
%
%  For the details, see:
%  M. Montardini, G. Sangalli, L. Tamellini, 
%  Optimal-order isogeometric collocation at Galerkin superconvergent points
%  Comput. Methods Appl. Mech. Engrg., 2016.
%
% Copyright (C) 2016, 2017 Monica Montardini
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

function col_pts = collocation_csp (knots, p)

if (p<3 || p>7)
     error ('The degree must be between 3 and 7')  
end 

% Superconvergent coordinates in parametric interval
jj7 = 0.5049185675126533;
jj5 = (sqrt(225-30*sqrt(30)))/15;
jj3 = 1/sqrt(3);

 flag = 0;
  if (~iscell (knots))
    flag = 1;
    knots = {knots};
  end
   ndir = numel (knots);
  for idir = 1:ndir
    uniqueknots  = unique (knots{idir});
    du           = diff (uniqueknots);
    
    nel = length (uniqueknots)-1; %number of elements
    midpoints=uniqueknots(1:end-1)+du./2;
    
    three_col1 = midpoints-du.*jj3./2;
    three_col2 = midpoints+du.*jj3./2;
    five_col1 = midpoints-du.*jj5./2;
    five_col2 = midpoints+du.*jj5./2;
    seven_col1 = midpoints-du.*jj7./2;
    seven_col2 = midpoints+du.*jj7./2;
    
    col_pts{idir} = zeros(1,nel+p-2);
    if rem(nel,2)~=0
        if (p==3)
            col_pts{idir}(1:2:end) = three_col1(1:2:end);
            col_pts{idir}(2:2:end) = three_col2(1:2:end);
        elseif (p==4)   
            v = zeros(1,nel);
            v(1:2:end) = midpoints(1:2:end);
            v(2:2:end-1) = uniqueknots(2:2:end-2);
            col_pts{idir} = [v(1:2),midpoints(2),v(3:end-1),uniqueknots(end-1),v(end)];
        elseif (p==5)
            v = zeros(1,nel+1);
            v(1:2:end-1) = five_col1(1:2:end);
            v(2:2:end) = five_col2(1:2:end);
            col_pts{idir} = [v(1:2),five_col1(2),v(3:end-2),five_col2(end-1),v(end-1:end)];   
        elseif (p==6)
            v = zeros(1,nel);
            v(1:2:end) = midpoints(1:2:end);
            v(2:2:end-1) = uniqueknots(2:2:end-2);
            col_pts{idir} = [v(1:2),midpoints(2),uniqueknots(3),v(3:end-1),midpoints(end-1),uniqueknots(end-1),v(end)];
        else
            v = zeros(1,nel+1);
            v(1:2:end-1) = seven_col1(1:2:end);
            v(2:2:end) = seven_col2(1:2:end);
            col_pts{idir} = [v(1:2),seven_col1(2),seven_col2(2),v(3:end-2),seven_col1(end-1),seven_col2(end-1),v(end-1:end)];
        end
    else
        if (p==3)
            col_pts{idir}(1:2:end-2) = three_col1(1:2:end-1);
            col_pts{idir}(2:2:end-1) = three_col2(1:2:end-1);
            col_pts{idir}(end) = three_col2(end);
        elseif (p==4)   
            v = zeros(1,nel);
            v(1:2:end) = midpoints(1:2:end-1);
            v(2:2:end) = uniqueknots(2:2:end-1);
            col_pts{idir} = [v(1:2),midpoints(2),v(3:end),midpoints(end)];
        elseif (p==5)
            v = zeros(1,nel);
            v(1:2:end-1) = five_col1(1:2:end-1);
            v(2:2:end) = five_col2(1:2:end-1);
            col_pts{idir} = [v(1:2),five_col1(2),v(3:end),five_col1(end),five_col2(end)];   
        elseif (p==6)
            v = zeros(1,nel);
            v(1:2:end) = midpoints(1:2:end-1);
            v(2:2:end) = uniqueknots(2:2:end-1);
            col_pts{idir}=[v(1:2),midpoints(2),uniqueknots(3),v(3:end-2),uniqueknots(end-2),v(end-1:end),midpoints(end)];
        else
            v = zeros(1,nel);
            v(1:2:end-1) = seven_col1(1:2:end-1);
            v(2:2:end) = seven_col2(1:2:end-1);
            col_pts{idir} = [v(1:2),seven_col1(2),seven_col2(2),v(3:end-2),seven_col2(end-2),v(end-1:end),seven_col1(end),seven_col2(end)];
        end
    end
  end
       
end