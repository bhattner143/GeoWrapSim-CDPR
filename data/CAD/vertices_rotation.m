function [vertices] = vertices_rotation(x_axis, y_axis, z_axis, input_vertices)
    rot = eul2rotm([x_axis, y_axis, z_axis],'XYZ') ;   % overall rotation matrix
    vertices = [rot*input_vertices']';

end