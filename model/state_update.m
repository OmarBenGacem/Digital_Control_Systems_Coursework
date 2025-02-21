function x_dot = state_update(x, u)

    s_dot  = x(2);
    s_ddot = - (F / M) * x(1) + u ;
    

    phi_dot  = x(4);
    phi_ddot = ( g/L + (F/(m * L) ) * x(1) - u/L ) * sin(x(3));

    x_dot = [s_dot ; s_ddot ; phi_dot, phi_ddot];

end