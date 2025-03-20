%% B4
%  Suppose that the control law is implemented with a discrete-time controller connected to the linearized
%  system (2) via an impulsive sampler (sampling the continuous-time state x(t) on one side) and a hold
%  (generating the continuous-time input u(t) on the other side). Find a value of the sampling time T
%  for which the closed-loop system is asymptotically stable. Find a value of T for which the closed-loop
%  system is unstable. Plot the corresponding y(t) and comment on the plots.


load_variables;
part_a;





%% Variables
t = sym("t");
k = sym('k');
K = sym("K", [1,4]);
Ts = sym("Ts");
tau = sym("tau");
kd = sym('Kd');
Kd = sym("Kd", [1,4]);

%% Functions
xd = symfun(sym("xd"), k); % digital representation of x
x = symfun(sym("x", [4,1]), t);
% dx = diff(x(t), t);
dx = symfun(sym("dx", [4,1]), t);
u = symfun(sym("u"), t);
ud = symfun(sym("ud"), k);




%% Code

lin_sys = dx == A * x + B * u;

cntrld_sys = dx == A * x + B * K * x;


%% DELETE - Is Repeated
% A_val_sym = subs(A, [F, g, M, L], [F_val, g_val, M_val, L_val]);
% B_val_sym = subs(B, [M, L], [M_val, L_val]);

% A_cl = subs(A, [F, g, M, L], [F_val, g_val, M_val, L_val]) + subs(B, [M, L], [M_val, L_val]) * subs(K, K, K_val)

T = sym("T");
Ad = sym("Ad", size(A));
Bd = sym("Bd", size(B));

F = A - B * Kd;
s = sym("s", [1, 4]);
lap = laplace(F);
eigz = diag(s) - F;
inverse = inv(eigz);
invlap = ilaplace(inverse)

