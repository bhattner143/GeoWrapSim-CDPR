function [t, int_flag] = ObtainCableObjectInterference(object_prop,A,P)
%UNTITLED12 Summary of this function goes here
%   Detailed explanation goes here
%Determine intersection PA and cylinders
    PA = A - P;
    syms t
    xx = (P(1) + t*PA(1));
    yy = (P(2) + t*PA(2));
    zz = (P(3) + t*PA(3));

    eqn = (xx - object_prop.center_base(1)).^2 + (zz - object_prop.center_base(3)).^2 - object_prop.a.^2==0;
    t = double(solve(eqn,t,'Real',true));
    
    if isempty(t)
        int_flag = false;
        % int_flag_array(index)          = 0;
    else
        int_flag = true;
        % int_flag_array(index)          = 1;
    end
end