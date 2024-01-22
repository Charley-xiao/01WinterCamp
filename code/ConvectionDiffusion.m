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
method = 1; % 1: central difference; 2: backward in convection;

Pe = rho * u * L / Tau;

% initial mesh
dx = L / num_nodes;
x = 0:dx:L;

% theoretical solution
phi_theory = phi_0 + (exp(x*Pe/L)-1)/(exp(Pe)-1)*(phi_L - phi_0);

% calculate the matrix
[A_E,A_W,A_P] = get_A(1,rho,u,dx,Tau);
[A_E2,A_W2,A_P2] = get_A(2,rho,u,dx,Tau);

% boundary condition
Q_0 = -A_E*phi_0;
Q_N = -A_E*phi_L;
phi(1)=phi_0;
phi2(1)=phi_0;
phi(num_nodes + 1)=phi_L;
phi2(num_nodes + 1)=phi_L;

% sovler linear equations
A=full(gallery('tridiag',num_nodes - 1,A_W,A_P,A_E));
A2=full(gallery('tridiag',num_nodes - 1,A_W2,A_P2,A_E2));
Q(1,1)=(Q_0);
Q(num_nodes - 1,1)=(Q_N);
Q(2:num_nodes - 2,1)=0;
phi(2:num_nodes)=A\Q;
phi2(2:num_nodes)=A2\Q;

plot(x,phi,'--or','LineWidth',2);
hold on;
plot(x,phi_theory,'-k','LineWidth',2);
hold on;
plot(x,phi2,'-x','LineWidth',2);
hold off
legend('central', 'theoretical', 'backward')
xlabel('x')
ylabel('\Phi')

function [A_E,A_W,A_P] = get_A(method,rho,u,dx,Tau)
    if (method == 2)
        A_E = min(rho*u,0)/dx-Tau/dx/dx;
        A_W = -max(rho*u,0)/dx-Tau/dx/dx;
        A_P = -A_E-A_W;
    elseif (method == 1)
        A_E = rho*u/dx-2*Tau/dx/dx;
        A_W = -rho*u/dx-2*Tau/dx/dx;
        A_P = -A_E-A_W;
    end
end

