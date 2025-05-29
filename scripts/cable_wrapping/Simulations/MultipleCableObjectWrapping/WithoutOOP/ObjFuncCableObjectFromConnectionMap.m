function f = ObjFuncCableObjectFromConnectionMap(bk_obs,object_connection_map)
    %UNTITLED14 Summary of this function goes here
    % The cable pattern used here 
    % A--straight--C3<--cwg3<--D3--st--C2<--cwg2<--D2--st--C1<--cwg1<--D1--st--P 
    % where < represents cable start to end direction

    num_of_obj_in_map = length(object_connection_map);%one less since one object is P
    
    [alpha_cell, C_cell, D_cell, delta_alpha_t_start_unit_cell, delta_alpha_t_end_unit_cell] =...
        GenerateCWGfromConnectionMap(object_connection_map, bk_obs);
    
    %Generate unit vectors along the straight part of cables and determine
    %the dot product with the tangent vector at the cable end points
    
    % For first object
    P = object_connection_map(1).object.P;% first object, point P
    A = object_connection_map(end).object.A;% first object, point P
    
    % Generated bewtween first object P and second object in object_connection_map
    DP_unit = (P - D_cell{object_connection_map(2).object.number})/norm(P - D_cell{object_connection_map(2).object.number}); %first (pt P) to last object
    f_start = norm(delta_alpha_t_start_unit_cell{object_connection_map(2).object.number}'* DP_unit - 1); 
    
    % For last object in object_connection_map
    CA_unit = (A   - C_cell{object_connection_map(end-1).object.number})/...
        norm(A    - C_cell{object_connection_map(end-1).object.number}); % last object to A        
    f_end   = norm(delta_alpha_t_end_unit_cell{object_connection_map(end-1).object.number}'* CA_unit - 1);
    
    % For middle objects
    f_middle = 0;

    for index = 2:num_of_obj_in_map-2 %num_of_obj_in_map-1 % Not required for first and last object, hence -2

        C_this_obj = C_cell{object_connection_map(index).object.number};
        for index2 = index+1:num_of_obj_in_map
            if isempty(D_cell{object_connection_map(index2).object.number}) == false
                D_next_obj = D_cell{object_connection_map(index2).object.number};
                C1D2_unit = (D_next_obj  - C_this_obj)/norm(D_next_obj - C_this_obj);
                D2C1_unit = -C1D2_unit;
        
                % dot product of current object cable end and vector C1D2
                f_middle  = f_middle + norm(delta_alpha_t_end_unit_cell{object_connection_map(index).object.number}'*C1D2_unit - 1);
                % dot product of next object cable start and vector D2C1
                f_middle  = f_middle + norm(delta_alpha_t_start_unit_cell{object_connection_map(index2).object.number}'*(D2C1_unit) - 1);

                break
            end
        end
    end
    f = f_start + f_end+ f_middle;  
end