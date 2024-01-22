clear all
clc

rho = 1;
L = 1;
u = 1;
Tau = 0.02;
phi_0 = 0;
phi_L = 1;
num_nodes = 50;
dt = 0.02;
max_t = 2;

dx = L / num_nodes;
x = 0:dx:L;

Pe = rho * u * L / Tau;
phi_theory = phi_0 + (exp(x*Pe/L)-1)/(exp(Pe)-1)*(phi_L - phi_0);

D = Tau * dt / (rho * dx^2);
C = u * dt / (2 * dx);

phi = linspace(phi_0, phi_L, num_nodes + 1)';

A_E = C - D;
A_P = 1 + 2 * D;
A_W = -C - D;

Q_0 = -A_W * phi_0;
Q_N = -A_E * phi_L;

t = 0.00;
figure;
while t < max_t
    phi(1) = phi_0;
    phi(end) = phi_L;
    
    A = diag(A_P * ones(num_nodes + 1, 1)) + diag(A_E * ones(num_nodes, 1), 1) + diag(A_W * ones(num_nodes, 1), -1);
    
    A(1, 1) = 1;
    A(end, end) = 1;

    b = phi;
    b(1,1) = (Q_0);
    b(end,1) = (Q_N);
    
    phi = A \ b;
    
    t = t + dt;
    if t > 0
        plot(x,phi_theory,'-k','LineWidth',2, 'DisplayName', 'Analytical Solution');
        hold on;
        plot(x, phi, '--', 'LineWidth', 2, 'DisplayName', ['t = ' num2str(t)]);
        hold on;
        xlabel('x');
        ylabel('\Phi');
        title('1D Convection-Diffusion Equation - Implicit Scheme');
        legend('show');
        drawnow;
        pause(0.01);
        hold off;
    end
end

