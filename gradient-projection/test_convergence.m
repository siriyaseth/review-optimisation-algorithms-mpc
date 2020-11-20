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
tau = 1e-12; gamma = 1;
show = true; % show iteration information or not


%% Test number of iterations before satisfying threshold with varying alpha

x = [0.2;0;-0.4;0];

niters = [];
alphas = 0.01:0.01:1;
for alpha = alphas
    [P,q] = mpc2qp_compact(x,Ad,Bd,Q,R,N);
    [U_iters,k] = solve_qp_grad_proj_box(P, q, kron(ones(N,1),umin), kron(ones(N,1),umax), alpha, 1, tau, true, u0);
    if (k > 10000)
        % break when current alpha takes too many iterations without
        % satisfying threshold
        break;
    end
    niters = [niters k];
end
scatter(alphas(1:length(niters)),niters,'LineWidth',1.5)
xlabel('\alpha')
ylabel('Number of iterations')
