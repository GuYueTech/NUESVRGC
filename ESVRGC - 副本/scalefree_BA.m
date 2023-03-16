% Copyright (C) 2006 Ginestra Bianconi <gbiancon@ictp.trieste.it>
% Copyright (C) 2007 Nicola Soranzo <soranzo@sissa.it>
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

% Create a directed scale-free network of n nodes with average degree 
% avg_degree and indegree and outdegree distributions P(k) ~ k^(-gamma) for 
% some gamma in [2, +Inf) using the Barabasi-Albert algorithm with common 
% preferential attachment for edge heads and tails
function [Y] = scalefree_BA(n, avg_degree)
if nargin < 2
    error('2 arguments required')
end
if ~isscalar(n) || n < 0
	error('n must be a scalar >= 0')
end
if ~isscalar(avg_degree) || avg_degree < 0
	error('avg_degree must be a scalar >= 0')
end
Y = zeros(n);

% Generate a first random edge
x = ceil(n * rand(1));
edge_found = false;
while ~edge_found
	y = ceil(n * rand(1));
	if x ~= y
		edge_found = true;
	end
end
Y(x, y) = 1;
% prob will contain tails and heads of all edges, one time for each edge
prob = [x y];

for i = 2:round(n * avg_degree)
	filtered_prob = [];
	% select a random node, 50% a new tail, 50% a new head
	if mod(i, 2)
		% choose a random tail x and filter out of prob x and all the direct successor of x
		while isempty(filtered_prob)
			x = ceil(n * rand(1));
			filtered_prob = prob(prob ~= x);
			to_filter = find(Y(x, :));
			for j = 1:numel(to_filter)
				filtered_prob = filtered_prob(filtered_prob ~= to_filter(j));
			end
		end
	else
		% choose a random head y and filter out of prob y and all the direct predecessor of y
		while isempty(filtered_prob)
			y = ceil(n * rand(1));
			filtered_prob = prob(prob ~= y);
			to_filter = find(Y(:, y));
			for j = 1:numel(to_filter)
				filtered_prob = filtered_prob(filtered_prob ~= to_filter(j));
			end
		end
	end

	% select the opposite end of the new edge using preferential attachment,
	% i.e. the probability of choosing a node is proportional to its degree
	if mod(i, 2)
    	y = filtered_prob(ceil(numel(filtered_prob) * rand(1)));
    else
    	x = filtered_prob(ceil(numel(filtered_prob) * rand(1)));
    end

    % Add the new edge and update prob
    Y(x, y) = 1;
    prob = [prob x y];
end
