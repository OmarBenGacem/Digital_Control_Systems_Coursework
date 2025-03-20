%% Define symbolic variables
% syms s s_dot phi phi_dot      
% syms u F M m g L      

%% Varaibles
s = sym("s");
t = sym("t");
s_dot = sym("s_dot");
phi = sym("phi");
phi_dot = sym("phi_dot");
F = sym("F");
M = sym("M");
m = sym("m");
g = sym("g");
L = sym("L");
A = sym("A", [4,4]);

%% Functions
x = symfun(sym("x", [4,1]), t);
x_dot = diff(x(t), t);
u = symfun(sym("u"), t);


x = [s; s_dot; phi; phi_dot];


x_eq = [0; 0; 0; 0];
u_eq = 0; 


s_dot_eq   = s_dot;
s_ddot_eq  = - (F/M)*s_dot + u/M;
phi_dot_eq = phi_dot;
phi_ddot_eq= ( (g/L) * sin(phi) ) + ( (F/(M*L)) * s_dot * cos(phi) ) - ( u/(L*M) * cos(phi) ) ;

x_dot = [s_dot_eq; s_ddot_eq; phi_dot_eq; phi_ddot_eq];

A_sym = jacobian(x_dot, x);

B_sym = jacobian(x_dot, u);
A = subs(A_sym, [s, s_dot, phi, phi_dot, u], [x_eq.' u_eq]);
B = subs(B_sym, [s, s_dot, phi, phi_dot, u], [x_eq.' u_eq]);


C = [1 0 0 0;
     0 0 1 0];
D = zeros(2, 1);

save("a3", "A", "B", "C", "D");

R = sym("R");
O = sym("O");
Reachability_Matrix = [B , A*B , A^2*B , A^3*B];
Observability_Matrix = [C ; C*A ; C*A^2 ; C*A^3];
Linear_Dynamics = [A , B ; C , D];



toOverleaf(Reachability_Matrix, "reachability", "R");
toOverleaf(Observability_Matrix, "observability", "O");
toOverleaf(Linear_Dynamics, "linear_dynamics");
toOverleaf(det(Reachability_Matrix),"reachability_det", "det(R)")
% toOverleaf(det(Observability_Matrix),"observability_det")
toOverleaf(A, "A");
toOverleaf(B, "B");
% toOverleaf(C, "C");
% toOverleaf(D, "D");


