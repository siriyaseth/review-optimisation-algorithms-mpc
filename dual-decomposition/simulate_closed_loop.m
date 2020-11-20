clear all
close all


%% Continuous time system

% Discrete-time system matrices

Phi = cell(2,2);
Phi{1,1} = 1; Phi{1,2} = 0;
Phi{2,1} = 0; Phi{2,2} = 1;

Gamma = cell(2,2);
Gamma{1,1} = 0.5; Gamma{1,2} = -0.5;
Gamma{2,1} = -0.5; Gamma{2,2} = 0.5;
% Cd = C;
nx = size(Gamma{1,1},1); nu = size(Gamma{1,1},2);

%% Setup constraints, reference and weights for MPC

% weights
Q1 = 10; Q2 = 10;
R1 = 1; R2 = 1;

% references
x1set = 19;
x2set = 81;

% constraitns
u1min = 0; u1max = 5;
u2min = 0; u2max = 5;

x1min = 0; x1max = 80;
x2min = 0; x2max = 80;

% horizon
N = 2;

%% Run simulation

gamma = 1; tau = 1e-4; show = true; visualise = false;
alpha = 0.45;
Ts = 0.5;
Tfinal = 20;
nsteps = round(Tfinal/Ts);
times = 0:Ts:Tfinal;

x1 = 80;
x2 = 20;
clf;
hold on
xs1 = zeros(nx,nsteps+1);
xs1(:,1) = x1;
xs2 = zeros(nx,nsteps+1);
xs2(:,1) = x2;

implementedU1 = zeros(nu,nsteps);
implementedU2 = zeros(nu,nsteps);

for i = 1:nsteps
    [P1,P2,q1,q2,A1,A2,b,C1,d1,C2,d2] = mpc2dualdecomp(Phi,Gamma,Q1,Q2,R1,R2,x1set,x2set,u1min,u1max,u2min,u2max,x1min,x1max,x2min,x2max,N,x1,x2);
    [y1,y2] = dual_decomp(P1,P2,q1,q2,A1,A2,b,C1,d1,C2,d2,alpha,gamma,tau,show,visualise);
    
    U1 = y1((N+1)*nx+1:(N+1)*nx+1+nu-1,1); U2 = y2((N+1)*nx+1:(N+1)*nx+1+nu-1,1);
    x1 = Phi{1,1}*x1 + Gamma{1,1}*U1 + Gamma{1,2}*U2; % agent 1
    x2 = Phi{2,2}*x2 + Gamma{2,1}*U1 + Gamma{2,2}*U2; % agent 2

    implementedU1(:,i) = U1;
    implementedU2(:,i) = U2;
    xs1(:,i+1) = x1;
    xs2(:,i+1) = x2;

    i
end


%% Plot

close all
figure(1)
plot(times,xs1(:),'b','LineWidth',1.5)
hold on
plot(times,xs2(:),'r','LineWidth',1.5)
xlabel('Time (s)')
ylabel('Number of Chickens')
legend('x_1','x_2')

figure(2)
plot(times(1:end-1),implementedU1(:),'b','LineWidth',1.5)
hold on
plot(times(1:end-1),implementedU2(:),'r','LineWidth',1.5)
xlabel('Time (s)')
ylabel('Chicken Flow Rate')
legend('u_1','u_2')

figure(3)
plot(times,xs1(:) + xs2(:))
