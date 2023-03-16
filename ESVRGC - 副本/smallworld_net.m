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

% Create a directed small-world network of n nodes with outdegree k using
% p as the probability of link rewiring
function Y = smallworld_net(n, k, p)
if nargin < 3
    error('3 arguments required')
end
if ~isscalar(n) || n < 0
	error('n must be a scalar >= 0')
end
if ~isscalar(k) || k < 0
	error('k must be a scalar >= 0')
end
if ~isscalar(p) || ~isreal(p) || p < 0 || p > 1
	error('p must be a scalar real >= 0 and <= 1')
end

Y = zeros(n);
% Connect each node to its k nearest neighbours
for i = 1:n
    for j = i-floor(k/2):i+ceil(k/2)
        if j ~= i
        	j = mod(j-1, n) + 1; % ensure j is in 1:n
            Y(i, j) = 1;
        end
    end
end

for i = 1:n
    for j = i-floor(k/2):i+ceil(k/2)
        if j ~= i
        	j = mod(j-1, n) + 1; % ensure j is in 1:n
			if (rand(1) <= p)
				% find a new destination node
				new_j = j;
				while new_j == i || Y(i, new_j) == 1
					new_j = ceil(rand(1) * n);
				end
				% drop link (i, j) and add (i, new_j)
				Y(i, j) = 0;
				Y(i, new_j) = 1;
			end
		end
	end
end
