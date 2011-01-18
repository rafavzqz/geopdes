% MP_EXTRACT_SUBDOMAINS: create a multipatch geometry extracting a
% set of subdomains from a gigen multipatch geometry.
%
% [new_geometry, new_interfaces, new_boundaries] = mp_extract_subdomains (geometry, interfaces, boundaries, subdomains)
%
% INPUT :
%
%   geometry, interfaces, boundaries: initial multipatch geometry
%   subdomains: list of subdomains to etract
%
% OUTPUT:
%
%   new_geometry, new_interfaces, new_boundaries: extracted multipatch geometry
%
% Copyright (C) 2010 Carlo de Falco, Rafael Vazquez
% 
% This program is free software; you can redistribute it and/or modify
% it under the terms of the GNU General Public License as published by
% the Free Software Foundation; either version 2 of the License, or
% (at your option) any later version.

% This program is distributed in the hope that it will be useful,
% but WITHOUT ANY WARRANTY; without even the implied warranty of
% MERCHANTABILITY or FITNESS FOR A PARTICULAR PURPOSE.  See the
% GNU General Public License for more details.
% 
% You should have received a copy of the GNU General Public License
% along with Octave; see the file COPYING.  If not, see
% <http://www.gnu.org/licenses/>.


function [new_geometry, new_interfaces, new_boundaries, new_cpt] = mp_extract_subdomains (geometry, interfaces, boundaries, subdomains)

npatch = numel (geometry);
nsubdomains = numel (subdomains);

new_subdomains = 1:npatch;
new_subdomains(subdomains) = 1:nsubdomains;

new_geometry   = geometry(subdomains);

keep_int = ismember ([interfaces.patch1], subdomains) & ismember ([interfaces.patch2], subdomains);
new_num_1 = num2cell (new_subdomains([interfaces.patch1]));
new_num_2 = num2cell (new_subdomains([interfaces.patch2]));
new_interfaces = interfaces;
[new_interfaces.patch1] = deal (new_num_1{:});
[new_interfaces.patch2] = deal (new_num_2{:});
new_interfaces = new_interfaces(keep_int);

keep_bnd = [];
nbnd     = numel (boundaries);
new_boundaries = boundaries;
for ibnd = 1:nbnd
  keep_patches = ismember (boundaries(ibnd).patches, subdomains);
  if (any (keep_patches));
    keep_bnd = [keep_bnd, ibnd];
    new_boundaries(ibnd).patches = new_subdomains (new_boundaries(ibnd).patches);
    new_boundaries(ibnd).patches = new_boundaries(ibnd).patches(keep_patches);
    new_boundaries(ibnd).nsides  = numel (new_boundaries(ibnd).patches);
    new_boundaries(ibnd).faces   = new_boundaries(ibnd).faces(keep_patches);
    new_boundaries(ibnd).flag    = ones (1, new_boundaries(ibnd).nsides);
    new_boundaries(ibnd).ornt1   = ones (1, new_boundaries(ibnd).nsides);
    new_boundaries(ibnd).ornt2   = ones (1, new_boundaries(ibnd).nsides);
  end
end
%new_boundaries = new_boundaries(keep_bnd);
%nbnd   = numel (new_boundaries);

int_2_bnd_1 = ismember ([interfaces.patch1], subdomains) & ~ismember ([interfaces.patch2], subdomains);
naddbnd = sum (int_2_bnd_1);
if (naddbnd)
  [new_boundaries(nbnd+(1:naddbnd)).name]     = deal (interfaces(int_2_bnd_1).ref);
  [new_boundaries(nbnd+(1:naddbnd)).nsides]   = deal (1);
  new_num  = num2cell (new_subdomains([interfaces(int_2_bnd_1).patch1]));
  [new_boundaries(nbnd+(1:naddbnd)).patches]  = deal (new_num{:});
  [new_boundaries(nbnd+(1:naddbnd)).faces]    = deal (interfaces(int_2_bnd_1).side1);
  for iii = nbnd+(1:naddbnd)
    [new_boundaries(iii).flag]     = 1;
    [new_boundaries(iii).ornt1]    = 1;
    [new_boundaries(iii).ornt2]    = 1;
  end
end

nbnd   = numel (new_boundaries);
int_2_bnd_2 = ~ismember ([interfaces.patch1], subdomains) & ...
    ismember ([interfaces.patch2], subdomains);
naddbnd = sum (int_2_bnd_2);
if (naddbnd)
  [new_boundaries(nbnd+(1:naddbnd)).name]     = deal (interfaces(int_2_bnd_2).ref);
  [new_boundaries(nbnd+(1:naddbnd)).nsides]   = deal (1);
  new_num  = num2cell (new_subdomains([interfaces(int_2_bnd_2).patch2]));
  [new_boundaries(nbnd+(1:naddbnd)).patches]  = deal (new_num{:});
  [new_boundaries(nbnd+(1:naddbnd)).faces]    = deal (interfaces(int_2_bnd_2).side2);
  [new_boundaries(nbnd+(1:naddbnd)).flag]     = deal (interfaces(int_2_bnd_2).flag);
  [new_boundaries(nbnd+(1:naddbnd)).ornt1]    = deal (interfaces(int_2_bnd_2).ornt1);
  [new_boundaries(nbnd+(1:naddbnd)).ornt2]    = deal (interfaces(int_2_bnd_2).ornt2);
end
end