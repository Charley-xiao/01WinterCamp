clear all
clc
% setup parameters
rho = 1;
L = 1;
u = 1;
Tau = 0.02;
phi_0= 0;
phi_L = 1;
num_nodes = 50;
method = 2; % 1: central difference; 2: backward in convection;

Pe = rho * u * L / Tau;

% intial mesh
dx = L / num_nodes;
x = 0:dx:L;

% theoretical solution
phi_theory = phi_0 + (exp(x*Pe/L)-1)/(exp(Pe)-1)*(phi_L - phi_0);

% calculate the matrix
if (method == 1)
    A_E = -2*Tau/(2*dx*dx) + (rho * u)/(2*dx);
    A_W = -2*Tau/(2*dx*dx) - (rho * u)/(2*dx);
    A_P = -(A_E + A_W);
elseif (method == 2)
    A_E = -2*Tau/(2*dx*dx);
    A_W = -2*Tau/(2*dx*dx) - (rho * u)/(dx);
    A_P = -(A_E + A_W);
end
% bounary condition
Q_0 = -A_W * phi_0;
Q_N = -A_E * phi_L;
phi(1)=phi_0;
phi(num_nodes + 1)=phi_L;

% sovler linear euqations
A=full(gallery('tridiag',num_nodes - 1,A_W,A_P,A_E));
Q(1,1)=(Q_0);
Q(num_nodes - 1,1)=(Q_N);
Q(2:num_nodes - 2,1)=0;
phi(2:num_nodes)=inv(A)*Q;

plot(x,phi,'--or','LineWidth',2);
hold on;
plot(x,phi_theory,'-k','LineWidth',2);
hold off
legend('calculated', 'theoretical')
xlabel('x')
ylabel('\Phi')

