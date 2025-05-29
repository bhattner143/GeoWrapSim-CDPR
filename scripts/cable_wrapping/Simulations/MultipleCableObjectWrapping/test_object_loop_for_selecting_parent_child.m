% Script file for Operational Space CTC
clc; clear; close all;
CASPR_log.SetLoggingDetails(CASPRLogLevel.INFO);

%%
objects = struct();
% 0--P, 1--O1, 2--O2, 3--O4, 4--O4, %--A
tot_num_objs = 5;
for obj_num = 1:tot_num_objs
    if obj_num == 1 % Cable attachment point P
        objects(obj_num).object.name   = 'P';
        objects(obj_num).object.number = 0;
        objects(obj_num).object.parent = 'P';
        objects(obj_num).object.child  = 'O4';
        objects(obj_num).object.child_number  = [];

    elseif obj_num < tot_num_objs && obj_num > 1 % Other objects
        objects(obj_num).object.name = strcat('O',num2str(obj_num-1));
        objects(obj_num).object.number = obj_num-1;
        objects(obj_num).object.parent = [];
        objects(obj_num).object.child  = [];
        objects(obj_num).object.child_number  = [];

    else
        objects(obj_num).object.name = strcat('O',num2str(obj_num-1));% for last object
        objects(obj_num).object.number = obj_num-1;
        objects(obj_num).object.parent = [];
        objects(obj_num).object.child  = 'A';
        objects(obj_num).object.child_number  = 5;

    end

end
%% Initial CWG is P-->04-->A
cwg_prev = 'PO4A';
% Initial object O1
obj_num_present = 1;

% Outer loop for looping through all objects
while obj_num_present < 4;
% for obj_num_present = 1:4
    %If the present object is a child of previous object then perform
    obj_num = 4;
    cwg_prev_inner_loop = cwg_prev;
    
    % Existing network map
    X = round([0 1 1 0])% Last three bits are the connection

    % Given an object indexed by obj_num_present, which object closest to
    % it, indexed by obj_num, is intersected by the CWG
    while obj_num > obj_num_present
        
        if X(obj_num) == 1
            intersection_flag = true;
        else
            intersection_flag = false;
        end
        % if there is a intersection then present cwg needs to be updated
        if intersection_flag
            %obj_num=4 points to 4th object, which is O3
            cwg_present = strcat(objects(obj_num).object.name,strcat('<--',cwg_prev_inner_loop));

            objects(obj_num_present).object.child         = objects(obj_num).object.name;
            objects(obj_num_present).object.child_number  = obj_num-1;
        
        % if there is no intersection then present cwg is same as
        % previous cwg
        else
            cwg_present = cwg_prev_inner_loop ;
%                 objects(obj_num_present).object.child = objects(5).object.name;
        end
    
        cwg_prev_inner_loop = cwg_present;
        obj_num             = obj_num - 1;
    end
    
    cwg_prev = strcat(cwg_prev ,strcat('<--',objects(obj_num_present).object.child))
    
    % Skip to the object which is connected to this object
    obj_num_present = objects(obj_num_present).object.child_number + 1;
    
    % For connecting the last object with O4
    if isempty(objects(obj_num_present).object.child)
        objects(obj_num_present).object.child = 'O4';
        objects(obj_num_present).object.child_number = 4;
    end
end


