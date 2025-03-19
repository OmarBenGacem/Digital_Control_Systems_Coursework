function x_dot = state_update(x, u)
    load_variables;

    
    s       = x(1);
    s_dot   = x(2);
    phi     = x(3);
    phi_dot = x(4);
    

    if false

        s_dot_eq   = max( s_dot, max_lin_v );
        s_ddot_eq  = - (F_val/M_val)*s_dot + u/M_val;
        phi_dot_eq = phi_dot;
        phi_ddot_eq= ( (g_val/L_val) * sin(phi) ) + ( (F_val/(M_val*L_val)) * s_dot * cos(phi) ) - ( u/(L_val*M_val) * cos(phi) ) ;
        
        
    else

        s_dot_eq   = s_dot;
        s_ddot_eq  = - (F_val/M_val)*s_dot + u/M_val;
        phi_dot_eq = phi_dot;
        phi_ddot_eq= ( (g_val/L_val) * sin(phi) ) + ( (F_val/(M_val*L_val)) * s_dot * cos(phi) ) - ( u/(L_val*M_val) * cos(phi) ) ;
        
    end

    x_dot = double([s_dot_eq; s_ddot_eq; phi_dot_eq; phi_ddot_eq]);
        
    

end