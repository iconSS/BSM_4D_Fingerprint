fid = fopen('graph.txt','rt');
g = {};
while ~feof(fid)
    thisline = fgetl(fid);
    disp(thisline)
    line = textscan(thisline, '%s %s %f %f %f %f');

    % With Inf distance, there is no connection
    % So we do not record them
    if ~isinf(line{end})
      g = [g; line]; 
    end
end

% finding all the names for the first column and second column
names1 = g(:,1);
names2 = g(:,2);
names = [names1(:)', names2(:)'];
% find all the unique names
unique_names = unique([names{:}]);

s_index = cellfun(@(x) find(contains(unique_names, x)), names1)';
t_index = cellfun(@(x) find(contains(unique_names, x)), names2)';

distances = g(:,6);
distances = [distances{:}];
% using the inverse of distances as the weights for the edges
% weights = 1 ./ distances;
% LWidths = 10 * weights ./ max(weights);

% creating the graph
G = graph(s_index, t_index, distances);
% using the inverse of distances as the weights for the edges
weights = 1 ./ G.Edges.Weight;
LWidths = 1 * (weights/10).^2; %./ max(weights);
plot(G,'NodeLabel', unique_names, 'MarkerSize',  10, 'NodeColor', 'r', 'LineWidth', LWidths, ...
       'EdgeColor', 'b', 'NodeFontSize', 10, 'EdgeLabel', G.Edges.Weight,'EdgeFontSize',16);