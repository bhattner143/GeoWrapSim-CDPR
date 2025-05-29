function f = ObjFuncForCableWrappingCASPR(b,k, obj_cable_wrapping, cable_num)

    % This objective function will work for both cylinder and cone since
    % the params generated from obj_cable_wrapping works accordingly with
    % the surface profile condition.

    optPram   = struct();

    psi_B     = 2*pi; 
    a         = obj_cable_wrapping.cable_wrapping_param.a_c;
    m         = obj_cable_wrapping.cable_wrapping_param.m_c;
    
    delta_k   = 0.001;
    % gradient of helix at pt B
    alpha1_k = cos(k*psi_B)*(a + k*m*psi_B);
    alpha2_k =                   -b*k*psi_B;
    alpha3_k = sin(k*psi_B)*(a + k*m*psi_B);
     
    alpha1_k_1 = cos((k-delta_k)*psi_B)*(a + (k-delta_k)*m*psi_B);
    alpha2_k_1 =                   -b*(k-delta_k)*psi_B;
    alpha3_k_1 = sin((k-delta_k)*psi_B)*(a + (k-delta_k)*m*psi_B);

    alpha_k_1 = [alpha1_k_1, alpha2_k_1, alpha3_k_1]';
    
    delta_alpha = [alpha1_k-alpha1_k_1, alpha2_k - alpha2_k_1, alpha3_k - alpha3_k_1]';
    delta_alpha_unit = delta_alpha/norm(delta_alpha);
    
    % PB_hat
    P = obj_cable_wrapping.cable_info.P_c(1:3);
    B = [alpha1_k, alpha2_k, alpha3_k]';
    BP_unit = (P - B)/norm(P - B);
    
    optPram.alpha_k_1 = alpha_k_1;
    optPram.delta_alpha_unit = delta_alpha_unit;
    optPram.B = B;
    optPram.BP_unit = BP_unit;

    f = norm(delta_alpha_unit'*BP_unit - 1);

end