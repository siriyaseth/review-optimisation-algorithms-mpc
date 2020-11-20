clear all
close all

%% System information

% Continuous-time system
Phic = [-0.0285 -0.0014; -0.0371 -0.1476];
Gammac = [-0.0850 0.0238; 0.0802 0.4462];
 
% Discretize CT system
Ts = 1.0;

Phi = expm(Phic*Ts);
syms tau
Gamma = double(int(expm(Phic*tau),0,Ts)*Gammac);
% Cd = C;

%% Setup MPC weights, constraints, reference and horizon

N = 3;

% u1min = -7; u1max = 7;
umin = [-5; -5];
umax = [5; 5];
xmin = [0; -0.4];
xmax = [4; 0.4];
nx = 2; nu = 2;

xset = [3;0]; % reference state

Q = [10 0; 0 0];
R = eye(2);

%% Augmented Lagrangian parameters

show = false; % show information each iteratino
visualise = false; % visualise solver

mu0 = 4; 
gamma1 = 1; % Increase rate for mu
tau = 1e-4;
alpha = 0.06;
gamma2 = 1;

%% Run simulation
Tfinal = 20;
nsteps = round(Tfinal/Ts);
times = 0:Ts:Tfinal;

x = [0;0];
clf;
hold on
xs = zeros(nx,nsteps+1);
xs(:,1) = x;

implementedU = zeros(nu,nsteps);
niters = [];
for i = 1:nsteps
    [P,q,A,b,C,d] = mpc2qp_box(Phi,Gamma,Q,R,xset,umin,umax,xmin,xmax,N,x);
    [Uauglag, k] = solve_qp_auglag_ineq(P,q,A,b,C,d,mu0,gamma1,tau,alpha,gamma2,show,visualise);
    U = Uauglag((N+1)*nx+1:(N+1)*nx+nu,1);

    x = Phi*x + Gamma*U;
    
    implementedU(:,i) = U;
    xs(:,i+1) = x;
    niters = [niters k];

end

%% Plot
close all
figure(1)
subplot(2,1,1)
plot(times,xs(1,:),'b','LineWidth',1.5)
ylabel('T (C)')
hold on
subplot(2,1,2)
plot(times,xs(2,:),'r','LineWidth',1.5)
ylabel('C_A (kgmol/m^3)')
xlabel('Time (s)')

figure(2)
subplot(2,1,1)
stairs(times(1:end-1),implementedU(1,:),'b','LineWidth',1.5)
ylabel('C_{Ai} (kgmol/m^3)')
hold on
subplot(2,1,2)
stairs(times(1:end-1),implementedU(2,:),'r','LineWidth',1.5)
ylabel('T_c (C)')
xlabel('Time (s)')

figure
plot(times(1:end-1),niters,'b','LineWidth',1.5);
xlabel('Time (s)')
ylabel('Number of Iterations')