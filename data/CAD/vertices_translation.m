function [vertices] = vertices_translation(T, input_vertices)

% This function translates every vertics in an imported STL object so the
% whole object is translated.

T = repmat(T, size(input_vertices, 1), 1);
vertices = input_vertices + T;

end