clear all
clc
% setup parameters
rho = 1;
L = 1;
u = 1;
Tau = 0.1;
phi_0= 0;
phi_L = 1;
num_nodes = 40;
max_t = 0.3;
dt = 0.003; %0.00325
method = 1; % 1: explicit Euler; 2: implicit Euler;

% intial mesh
num_t = max_t / dt;
dx = L / num_nodes;
x = 0:dx:L;
Q(1:num_nodes+1,1)=0;
phi(1:num_nodes,1)=0;
phi(num_nodes + 1) = phi_L;

ih = 1;
window_ouput = cast(num_t/5,"uint8");
for it = 1:num_t
    % bounary condition
    phi(1)=phi_0;
    phi(num_nodes + 1)=phi_L;
    
    if (method == 1)
        A_E = Tau/rho/(dx*dx)*dt - u/(2*dx)*dt;
        A_W = Tau/rho/(dx*dx)*dt + u/(2*dx)*dt;
        A_P = 1 - 2*Tau/rho/(dx*dx)*dt;
        Q = phi;
        A=full(gallery('tridiag',num_nodes + 1,A_W,A_P,A_E));
        phi=A*Q;
    elseif (method == 2)
        A_E = -Tau/rho/(dx*dx)*dt + u/(2*dx)*dt;
        A_W = -Tau/rho/(dx*dx)*dt - u/(2*dx)*dt;
        A_P = 1 + 2*Tau/rho/(dx*dx)*dt;
        Q = phi;
        A=full(gallery('tridiag',num_nodes + 1,A_W,A_P,A_E));
        phi=A\Q;
    end

    if mod(it, window_ouput) == 0
        phi(1)=phi_0;
        phi(num_nodes + 1)=phi_L;
        h(ih) = plot(x,phi,'--o','LineWidth',2,...
            'DisplayName', sprintf('t =  %.2f', it*dt));
        ih=ih+1;
        hold on;
    end
end

hold off
hL = legend(h);
set(hL,'Location','west')
xlabel('x')
ylabel('\Phi')

