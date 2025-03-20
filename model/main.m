close all;
clc;
clear;
addpath("analytic_work\");
part_a;

load_variables;

plot_B2  = true;
plot_B3  = true;
plot_B4  = true;
plot_B5  = true;
plot_B7  = true;
plot_B8  = true;
plot_B9  = true;
plot_B10 = true;

%                      s     sdot   phi     dphi
initial_conditions = [ 0,   0,     0.0872665,      0;
                       0,     0.1,   -0.174533, 0;   
                       0,     0,     -0.261799,   0;  
                       0,     0,     0.05,   0];  
tspan = [0 10];








%% B1 B2
%  Design a state feedback controller u(t) = Kx(t) which stabilizes the linearized system (2) in x(t) = 0. Comment your design



A_val_sym = subs(A, [F, g, M, L], [F_val, g_val, M_val, L_val]);
B_val_sym = subs(B, [M, L], [M_val, L_val]);

A_val = double(A_val_sym);
B_val = double(B_val_sym);

% poles = [-1, -190, -54, -20]
% poles = [-1, -2, -4, -8];
poles=[-0.5,-1,-1.5,-2];
% toOverleaf(P==poles, "poles");
K = place(A_val, B_val, poles);

A_cl = A_val - B_val*K;


toOverleaf(poles, "poles", true)
toOverleaf(K, "K", true)

if plot_B2

    figure('Position', [100, 100, 2200, 800]);
    nSim = size(initial_conditions, 1);
    
    ax1_all = [];  % Store subplot handles for linking
    ax2_all = [];
    
    for i = 1:nSim
        x0 = initial_conditions(i, :).';
        [t, x] = ode45(@(t, x) A_cl*x, tspan, x0);
    
        % Convert initial condition from radians to degrees
        IC_deg = initial_conditions(i,3) * (180/pi);
    
        % First subplot: State evolution
        ax1 = subplot(2, nSim, i);
        yyaxis left
        plot(t, x(:,1)*100, '-', 'LineWidth', 1.5)
        ylabel('s (cm)')
    
        % Find symmetric limits for 's (cm)'
        y_left_min = min(x(:,1)*100);
        y_left_max = max(x(:,1)*100);
        y_left_lim = max(abs([y_left_min, y_left_max]));  % Symmetric limit
        ylim([-y_left_lim, y_left_lim]);  % Ensure zero alignment
    
        yyaxis right
        plot(t, x(:,3)*180/pi, '--', 'LineWidth', 1.5)
        ylabel('\phi (deg)')
    
        % Find symmetric limits for 'φ (deg)'
        y_right_min = min(x(:,3)*180/pi);
        y_right_max = max(x(:,3)*180/pi);
        y_right_lim = max(abs([y_right_min, y_right_max]));  % Symmetric limit
        ylim([-y_right_lim, y_right_lim]);  % Ensure zero alignment
    
        title(sprintf("Simulation %d, (IC: %.2f°)", i, IC_deg))  % Updated title with degree symbol
        xlabel('Time (s)')
        grid on;
        legend('s (cm)', '\phi (deg)', 'Location', 'Best')
        xlim([0 t(end)]); 
        ax1_all = [ax1_all, ax1];  % Store handle
    
        % Second subplot: Actuation input
        ax2 = subplot(2, nSim, i + nSim);
        u = K * x';  % Ensure matrix multiplication is correct
        plot(t, u, 'LineWidth', 1.5)
        title(sprintf("Simulation %d Actuation, (IC: %.2f°)", i, IC_deg))  % Updated title with degree symbol
        xlabel('Time (s)')
        ylabel('u(t)')
        grid on;
        xlim([0 t(end)]);
        ax2_all = [ax2_all, ax2];  % Store handle
    end
    
    % Link x-axes for better visualization
    linkaxes([ax1_all, ax2_all], 'x');


    sgtitle('B2) Continuous Time System Responses: Linear and Angular Displacement')
    saveas(gcf, '../figures/b2.png');
end



%% B3
%  Display plots of y(t) for the nonlinear system (1), from the same initial states x(0) and using the controller designed in point B1. Comment on these plots.

if plot_B3


    figure('Position', [100, 100, 2200, 800]);
    nSim = size(initial_conditions, 1);
    
    ax1_all = [];  % Store subplot handles for linking
    ax2_all = [];
    
    for i = 1:nSim
        x0 = initial_conditions(i, :).';
        [t, x] = ode45(@(t, x) state_update(x, -K*x), tspan, x0);
    
        % Convert initial condition from radians to degrees
        IC_deg = initial_conditions(i,3) * (180/pi);
    
        % First subplot: State evolution
        ax1 = subplot(2, nSim, i);
        yyaxis left
        plot(t, x(:,1)*100, '-', 'LineWidth', 1.5)
        ylabel('s (cm)')
    
        % Find symmetric limits for 's (cm)'
        y_left_min = min(x(:,1)*100);
        y_left_max = max(x(:,1)*100);
        y_left_lim = max(abs([y_left_min, y_left_max]));  % Symmetric limit
        ylim([-y_left_lim, y_left_lim]);  % Ensure zero alignment
    
        yyaxis right
        plot(t, x(:,3)*180/pi, '--', 'LineWidth', 1.5)
        ylabel('\phi (deg)')
    
        % Find symmetric limits for 'φ (deg)'
        y_right_min = min(x(:,3)*180/pi);
        y_right_max = max(x(:,3)*180/pi);
        y_right_lim = max(abs([y_right_min, y_right_max]));  % Symmetric limit
        ylim([-y_right_lim, y_right_lim]);  % Ensure zero alignment
    
        title(sprintf("Simulation %d, (IC: %.2f°)", i, IC_deg))  % Updated title with degree symbol
        xlabel('Time (s)')
        grid on;
        legend('s (cm)', '\phi (deg)', 'Location', 'Best')
        xlim([0 t(end)]); 
        ax1_all = [ax1_all, ax1];  % Store handle
    
        % Second subplot: Actuation input
        ax2 = subplot(2, nSim, i + nSim);
        u = K * x';  % Ensure matrix multiplication is correct
        plot(t, u, 'LineWidth', 1.5)
        title(sprintf("Simulation %d Actuation, (IC: %.2f°)", i, IC_deg))  % Updated title with degree symbol
        xlabel('Time (s)')
        ylabel('u(t)')
        grid on;
        xlim([0 t(end)]);
        ax2_all = [ax2_all, ax2];  % Store handle
    end
    
    % Link x-axes for better visualization
    linkaxes([ax1_all, ax2_all], 'x');


    sgtitle('B3) Nonlinear System Responses: Displacement and Angle')
    saveas(gcf, '../figures/b3_x.png');



    figure('Position', [100, 100, 2200, 800]);
    nSim = size(initial_conditions, 1);
    
    ax1_all = [];  % Store subplot handles for linking
    ax2_all = [];
    
    for i = 1:nSim
        x0 = initial_conditions(i, :).';
        [t, x] = ode45(@(t, x) state_update(x, -K*x), tspan, x0);
    
        % Convert initial condition from radians to degrees
        IC_deg = initial_conditions(i,3) * (180/pi);
    
        % First subplot: State evolution
        ax1 = subplot(2, nSim, i);
        yyaxis left
        plot(t, x(:,2)*100, '-', 'LineWidth', 1.5)
        ylabel('s (cm)')
    
        % Find symmetric limits for 's (cm)'
        y_left_min = min(x(:,2)*100);
        y_left_max = max(x(:,2)*100);
        y_left_lim = max(abs([y_left_min, y_left_max]));  % Symmetric limit
        ylim([-y_left_lim, y_left_lim]);  % Ensure zero alignment
    
        yyaxis right
        plot(t, x(:,4)*180/pi, '--', 'LineWidth', 1.5)
        ylabel('\phi (deg)')
    
        % Find symmetric limits for 'φ (deg)'
        y_right_min = min(x(:,4)*180/pi);
        y_right_max = max(x(:,4)*180/pi);
        y_right_lim = max(abs([y_right_min, y_right_max]));  % Symmetric limit
        ylim([-y_right_lim, y_right_lim]);  % Ensure zero alignment
    
        title(sprintf("Simulation %d, (IC: %.2f°)", i, IC_deg))  % Updated title with degree symbol
        xlabel('Time (s)')
        grid on;
        legend('s (cm)', '\phi (deg)', 'Location', 'Best')
        xlim([0 t(end)]); 
        ax1_all = [ax1_all, ax1];  % Store handle
    
        % Second subplot: Actuation input
        ax2 = subplot(2, nSim, i + nSim);
        u = K * x';  % Ensure matrix multiplication is correct
        plot(t, u, 'LineWidth', 1.5)
        title(sprintf("Simulation %d Actuation, (IC: %.2f°)", i, IC_deg))  % Updated title with degree symbol
        xlabel('Time (s)')
        ylabel('u(t)')
        grid on;
        xlim([0 t(end)]);
        ax2_all = [ax2_all, ax2];  % Store handle
    end
    
    % Link x-axes for better visualization
    linkaxes([ax1_all, ax2_all], 'x');

    sgtitle('B3) Nonlinear System Responses: Angular and Linear Velocities')
    saveas(gcf, '../figures/b3_v.png');

end





%% B4
%  Suppose that the control law is implemented with a discrete-time controller connected to the linearized
%  system (2) via an impulsive sampler (sampling the continuous-time state x(t) on one side) and a hold
%  (generating the continuous-time input u(t) on the other side). Find a value of the sampling time T
%  for which the closed-loop system is asymptotically stable. Find a value of T for which the closed-loop
%  system is unstable. Plot the corresponding y(t) and comment on the plots.


% LECTURE 6
% LECTURE 7 SLIDE 17/19

part_b_4;
Ts = 0.2;




%% B5
%  Similarly to point B4, apply the discrete-time control law to the nonlinear system (1). This time,
%  focus and discuss about the possible degradation in performance for different values of T. Display
%  plots of y(t) and comment on the plots.

if plot_B5

end


%% B6
%  A bunch of equations


Ad = expm((A_val * Ts));
% syms tau
% Bd = int(expm(A_val * (Ts - tau)) * B_val, tau, [0, Ts]);
Bd = integral( @(tau) expm(A_val * (Ts - tau)  ) * B_val, 0, Ts, 'ArrayValued', true );
Cd = C;
Dd = D;

toOverleaf(Ad, "Ad", true);
toOverleaf(Bd, "Bd", true);
toOverleaf(Cd, "Cd", true)
toOverleaf(Dd, "Dd", true);


%% B7
%  Design a discrete-time state feedback control law u(k) = Kdx(k) such that the closed-loop associated
%  to the discretized linearized system computed in point B6 is asymptotically stable (Hint: place the
%  poles of Ad + KdBd mapping the eigenvalues you chose in point B1 to the z-plane). Display plots
%  of y(t) and comment on the plots.


Kd = sym("Kd", [1,4]);
poles_discrete = exp(poles * Ts);
Kd = place(Ad, Bd, poles_discrete);

A_d_cl = Ad - Bd * Kd; 


toOverleaf(poles_discrete, "poles_d", true);
toOverleaf(Kd, "Kd", true);



if plot_B7
    
    figure('Position', [100, 100, 1600, 800]);
    for i = 1:size(initial_conditions,1)


        state = initial_conditions(i, :).';  % Column vector

        t_sim = tspan(1):Ts:tspan(end);
        x_hist = zeros(4, length(t_sim)); 
        
        % Simulate the discrete-time system
        for idx = 1:length(t_sim)
            x_hist(:, idx) = state;
            state = A_d_cl*(state);
        end
        

        subplot(2,2,i)
        yyaxis left
        plot(t_sim, x_hist(1, :)*100, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4) % Discrete points
        ylabel('s (cm)')

        yyaxis right
        plot(t_sim, x_hist(3, :)*180/pi, 's--', 'LineWidth', 1.5, 'MarkerSize', 4) % Discrete points
        ylabel('\phi (deg)')
        
        title(sprintf('Initial Condition %d', i))
        xlabel('Time (s)')
        grid on;
        legend('s (cm)', '\phi (deg)', 'Location', 'Best')
    end
    
    sgtitle('Discrete-Time System Responses: Displacement and Angle')
    saveas(gcf, '../figures/Discrete_Time_System_responce_x.png');


    figure('Position', [100, 100, 1600, 800]);
    for i = 1:size(initial_conditions,1)


        state = initial_conditions(i, :).';
        

        t_sim = tspan(1):Ts:tspan(end);
        x_hist = zeros(4, length(t_sim));

        for idx = 1:length(t_sim)
            x_hist(:, idx) = state; 
            state = A_d_cl*(state); 
        end
        

        subplot(2,2,i)
        yyaxis left
        plot(t_sim, x_hist(2, :)*100, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4)
        ylabel('Linear Velocity (cm/s)')

        yyaxis right
        plot(t_sim, x_hist(4, :)*180/pi, 's--', 'LineWidth', 1.5, 'MarkerSize', 4)
        ylabel('Angular Velocity (deg/s)')
        
        title(sprintf('Initial Condition %d', i))
        xlabel('Time (s)')
        grid on;
        legend('cm/s', 'rad/s', 'Location', 'Best')
    end
    
    sgtitle('B7) Discrete-Time System Responses: Linear and Angular Velocities')
    saveas(gcf, '../figures/Discrete_Time_System_responce_v.png');



end




%% B8
%  Apply the discrete-time state feedback control law u(k) = Kdx(k) to the continuous-time nonlinear
%  system (1) (obviously sampling the state and holding the input). Display plots of y(t), comment on
%  the plots and compare with the results of point B5.

if plot_B8
    figure('Position', [100, 100, 2200, 800]);
    nSim = size(initial_conditions, 1);
    
    for i = 1:nSim



        x0 = initial_conditions(i, :).';
        t_total = [];
        x_total = [];
        t_u = []; 
        u_total = []; 
        
        t_offset = 0;
        T_final = tspan(end); 
        N_intervals = ceil(T_final / Ts);
        
        for k = 1:N_intervals

            u = -Kd * x0;
            [t_interval, x_interval] = ode45(@(t, x) state_update(x, u), [0 Ts], x0);
           
            t_interval = t_interval + t_offset;

            t_total = [t_total; t_interval];
            x_total = [x_total; x_interval];

            t_u = [t_u; t_offset; t_offset + Ts];
            u_total = [u_total; u; u];

            x0 = x_interval(end, :)';
        
            t_offset = t_offset + Ts;
        end
      
        ax1 = subplot(2, nSim, i);
        yyaxis left
        plot(t_total, x_total(:,1)*100, '-', 'LineWidth', 1.5)
        ylabel('s (cm)')
        yyaxis right
        plot(t_total, x_total(:,3)*180/pi, '--', 'LineWidth', 1.5)
        ylabel('\phi (deg)')
        title(sprintf('IC %d State', i))
        xlabel('Time (s)')
        grid on;
        legend('s (cm)', '\phi (deg)', 'Location', 'Best')
        xlim([0 T_final]); 
        ax2 = subplot(2, nSim, i+nSim);
        stairs(t_u, u_total, 'LineWidth', 1.5)
        title(sprintf('IC %d Actuation', i))
        xlabel('Time (s)')
        ylabel('u(t)')
        grid on;
        xlim([0 T_final]);

        linkaxes([ax1, ax2], 'x');
    end
    
    sgtitle('B8) Continuous-Time System with Discrete-Time Control and Actuation u(t)')
    saveas(gcf, '../figures/Continuous_Time_Discrete_Actuation.png');
end







%% B8
%  Apply the discrete-time state feedback control law u(k) = Kdx(k) to the continuous-time nonlinear
%  system (1) (obviously sampling the state and holding the input). Display plots of y(t), comment on
%  the plots and compare with the results of point B5.





%% B9
%  A bunch of equations


if plot_B9
    
    figure('Position', [100, 100, 1600, 800]);
    for i = 1:size(initial_conditions,1)
        

        state = initial_conditions(i, :).';
        init = initial_conditions(i, :).';
        save("digital_values.mat","Ad","Bd", "tspan", "Ts", "init");
        Kd_star = fminsearch(@optimization_problem, Kd);
        delete("digital_values.mat")
        A_d_cl_star = ((Ad - Bd * Kd_star)); 
        t_sim = tspan(1):Ts:tspan(end);
        x_hist = zeros(4, length(t_sim));
        
        for idx = 1:length(t_sim)
            x_hist(:, idx) = state; 
            state = A_d_cl_star*(state); 
        end
        

        subplot(2,2,i)
        yyaxis left
        plot(t_sim, x_hist(1, :)*100, 'o-', 'LineWidth', 1.5, 'MarkerSize', 4) 
        ylabel('s (cm)')

        yyaxis right
        plot(t_sim, x_hist(3, :)*180/pi, 's--', 'LineWidth', 1.5, 'MarkerSize', 4)
        ylabel('\phi (deg)')


        title(sprintf('Initial Condition %d', i))
        xlabel('Time (s)')
        grid on;
        legend('s (cm)', '\phi (deg)', 'Location', 'Best')
    end

       sgtitle('B9) Discrete-Time System Responses Using an Optimal Controller: s (cm) and \phi (deg)')
       saveas(gcf, '../figures/Discrete_Time_Discrete_Actuation_Optimal.png');
    
end

function cost = optimization_problem(k)
    persistent Ad Bd ic tspan Ts Q J

    if isempty(Ad)
        load("digital_values.mat")
        Q = diag([100, 1, 100, 1]);
        J = @(x,u) x'*Q*x + u'*u;
        Ad = Ad;
        Bd = Bd;
        ic = init;
        tspan = tspan;
        Ts = Ts;
    end


    state = ic;
    t_sim = tspan(1):Ts:tspan(end);
    x_hist = zeros(4, length(t_sim));

    A_d_cl = ((Ad - Bd * k)); 

    for idx = 1:length(t_sim)
        x_hist(:, idx) = state; 
        state = A_d_cl*(state); 
    end

    cost= sum(arrayfun(@(i) J(x_hist(:, i), k * x_hist(:, i)), 1:size(x_hist,2)));
end



%% B10
%  Apply the discrete-time optimal control law u(k) = K∗dx(k) to the continuous-time nonlinear system
%  (1). Display plots of y(t), comment on the plots and compare with the results of point B5 and B8.
%  (Hint: in point B5 you applied a controller designed for continuous-time as a discrete-time controller.
%  In point B8 you applied a controller designed for discrete-time as a discrete-time controller. In point
%  B10 you applied an optimal controller designed for discrete-time as a discrete-time controller. Thus,
%  the performances of the three points should be different)


if plot_B10

    figure('Position', [100, 100, 2200, 800]);
    nSim = size(initial_conditions, 1);
    
    for i = 1:nSim
        % Get the i-th initial condition (as a column vector)
        x0 = initial_conditions(i, :).';

        init = x0;
        save("digital_values.mat","Ad","Bd", "tspan", "Ts", "init");
        Kd_star = fminsearch(@optimization_problem, Kd);
        delete("digital_values.mat")
        
        t_total = [];
        x_total = [];
        t_u = [];
        u_total = [];
        
        t_offset = 0;
        T_final = tspan(end);
        N_intervals = ceil(T_final / Ts);
        
        for k = 1:N_intervals
            u = -Kd_star * x0;
            
            [t_interval, x_interval] = ode45(@(t, x) state_update(x, u), [0 Ts], x0);
            t_interval = t_interval + t_offset;

            t_total = [t_total; t_interval];
            x_total = [x_total; x_interval];

            t_u = [t_u; t_offset; t_offset + Ts];
            u_total = [u_total; u; u];

            x0 = x_interval(end, :)';
            

            t_offset = t_offset + Ts;
        end

        ax1 = subplot(2, nSim, i);
        yyaxis left
        plot(t_total, x_total(:,1)*100, '-', 'LineWidth', 1.5)
        ylabel('s (cm)')
        yyaxis right
        plot(t_total, x_total(:,3)*180/pi, '--', 'LineWidth', 1.5)
        ylabel('\phi (deg)')
        title(sprintf('IC %d State', i))
        xlabel('Time (s)')
        grid on;
        legend('s (cm)', '\phi (deg)', 'Location', 'Best')
        xlim([0 T_final]);
        

        ax2 = subplot(2, nSim, i+nSim);
        stairs(t_u, u_total, 'LineWidth', 1.5)
        title(sprintf('IC %d Actuation', i))
        xlabel('Time (s)')
        ylabel('u(t)')
        grid on;
        xlim([0 T_final]);
        

        linkaxes([ax1, ax2], 'x');
    end
    
    sgtitle('B10) Continuous-Time System with an Optimal Discrete-Time Control and Actuation u(t)')
    saveas(gcf, '../figures/Continuous_Time_Discrete_Actuation_Optimal.png');
end










