% Copyright (C) 2006 Ginestra Bianconi <gbiancon@ictp.trieste.it>
% Copyright (C) 2006,2007 Nicola Soranzo <soranzo@sissa.it>
%
% This program is free software; you can redistribute it and/or modify it 
% under the terms of the GNU General Public License as published by the Free 
% Software Foundation; either version 2 of the License, or (at your option) 
% any later version.
%
% This program is distributed in the hope that it will be useful, but WITHOUT
% ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or 
% FITNESS FOR A PARTICULAR PURPOSE. See the GNU General Public License 
% for more details.

% You should have received a copy of the GNU General Public License 
% along with this program; see the file COPYING. If not, write to the Free 
% Software Foundation, Inc., 51 Franklin Street, Fifth Floor, Boston, MA 
% 02110-1301, USA.

% Create a directed random network of n nodes with average node degree 
% avg_degree
function Y = random_net(n, avg_degree)
if nargin < 2
    error('2 arguments required')
end
if ~isscalar(n) || n < 0
	error('n must be a scalar >= 0')
end
if ~isscalar(avg_degree) || avg_degree < 0
	error('avg_degree must be a scalar >= 0')
end

R = rand(n);
link_prob = avg_degree / (n-1);
Y = zeros(n);
for i = 1:n
     for j = 1:n
         if i ~= j && R(i, j) < link_prob
         	Y(i, j) = 1;
         end
     end
end
