%% Closed-loop simulation of linear tilted table system

clear all
close all
%% Define System

% Continuous-time linear system
g = 9.81;

A = [0 1 0 0;
    0 0 0 0;
    0 0 0 1;
    0 0 0 0];

B = [0                       0;
    -g*(5.0/7.0)            0;
    0                       0;
    0            -g*(5.0/7.0)];

C = [1 0 0 0;
    0 0 1 0];

% Discrete-time system

Ts = 1.0/25;

Ad = expm(A*Ts);
syms tau
Bd = double(int(expm(A*tau),0,Ts)*B);
Cd = C;

nx = 4; % Number of states
nu = 2; % Number of inputs

%% Setup constraints and weights

umax = [0.06; 0.06];
umin = [-0.06; -0.06];

rho = 10;
Q = rho*eye(4);
R = eye(2);
N = 7;


%% Optimisation algorithm parameters
alpha = 0.1;
tau = 1e-12; gamma = 1;
show = true; % show iteration information or not


%% Run our gradient projection algorithm with no warmstart
Tfinal = 15;
nsteps = round(Tfinal/Ts);
times = 0:Ts:Tfinal;
warmstart = false;

x = [0.4;-0.4;-0.45;-0.15];

xs_gradproj = zeros(nx,nsteps+1);
xs_gradproj(:,1) = x;
implementedU_gradproj = zeros(nu,nsteps);
niters = [];
for i = 1:nsteps
    if warmstart && (i > 1)
        u_warm = U_iters(:,end);
    else
        u_warm = zeros(N*nu,1);
    end
    [P,q] = mpc2qp_compact(x,Ad,Bd,Q,R,N);
    [U_iters,k] = solve_qp_grad_proj_box(P, q, kron(ones(N,1),umin), kron(ones(N,1),umax), alpha, gamma, tau, true, u_warm);
    U = U_iters(1:nu,end);
    
    tspan = [0 Ts];
    [t,xlist] = ode45(@(t,x) tilted_table_ode(t,x,U), tspan, x);
    x = xlist(end,:)';

    
    implementedU_gradproj(:,i) = U(:,1);
    xs_gradproj(:,i+1) = x;
    
    niters = [niters k];
end

%% View iterations requierd for each optimisation problem
figure
plot(times(1:end-1),niters,'b','LineWidth',1.5);
hold on
%% Run our gradient projection algorithm with warmstart
Tfinal = 15;
nsteps = round(Tfinal/Ts);
times = 0:Ts:Tfinal;
warmstart = true;

x = [0.4;-0.4;-0.45;-0.15];

implementedU_gradproj = [];
xs_gradproj = zeros(nx,nsteps+1);
xs_gradproj(:,1) = x;
implementedU_gradproj = zeros(nu,nsteps);
niters = [];
for i = 1:nsteps
    if warmstart && (i > 1)
        u_warm = U_iters(:,end);
    else
        u_warm = zeros(N*nu,1);
    end
    [P,q] = mpc2qp_compact(x,Ad,Bd,Q,R,N);
    [U_iters,k] = solve_qp_grad_proj_box(P, q, kron(ones(N,1),umin), kron(ones(N,1),umax), 0.1, 1, tau, true, u_warm);
    U = U_iters(1:nu,end);
    
    tspan = [0 Ts];
    [t,xlist] = ode45(@(t,x) tilted_table_ode(t,x,U), tspan, x);
    x = xlist(end,:)';

    
    implementedU_gradproj(:,i) = U(:,1);
    xs_gradproj(:,i+1) = x;
    
    niters = [niters k];
end
%% View iterations requierd for each optimisation problem
plot(times(1:end-1),niters,'r','LineWidth',1.5);

xlabel('Time (s)')
ylabel('Number of Iterations')
legend('No warmstart','Warmstart')
