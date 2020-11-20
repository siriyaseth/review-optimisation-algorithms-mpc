clear all
close all


%% Continuous time system

A = [0 0;
     0 0];

B = [1  -1;
     -1  1];
 
% Discretize CT system

Ts = 1.0/2;

Ad = expm(A*Ts);
syms tau
Bd = double(int(expm(A*tau),0,Ts)*B);
% Cd = C;


%% MPC chicken coop setup

Amat = cell(2,2);
Amat{1,1} = 1; Amat{1,2} = 0;
Amat{2,1} = 0; Amat{2,2} = 1;

Bmat = cell(2,2);
Bmat{1,1} = 0.5; Bmat{1,2} = -0.5;
Bmat{2,1} = -0.5; Bmat{2,2} = 0.5;

nx = size(Bmat{1,1},1); nu = size(Bmat{1,1},2);

Q1 = 10; Q2 = 10;
R1 = 1; R2 = 1;
% 
% x1set = 80;
% x2set = 20;
% 
% x1set = 19;
% x2set = 81;

x1set = 19;
x2set = 81;

u1min = 0; u1max = 5;
u2min = 0; u2max = 5;

x1min = 0; x1max = 80;
x2min = 0; x2max = 80;

N = 2;

x10 = 40;
x20 = 40;


%% Run simulation no warmstart

alpha0 = 0.1; gamma = 1; tau = 1e-4; show = false; visualise = false;
alpha0 = 0.45;
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
niters = [];
for i = 1:nsteps
    [P1,P2,q1,q2,A1,A2,b,C1,d1,C2,d2] = mpc2dualdecomp(Amat,Bmat,Q1,Q2,R1,R2,x1set,x2set,u1min,u1max,u2min,u2max,x1min,x1max,x2min,x2max,N,x1,x2);
    [y1,y2,info] = dual_decomp(P1,P2,q1,q2,A1,A2,b,C1,d1,C2,d2,alpha0,gamma,tau,show,visualise);

%     P = blkdiag(P1,P2);
%     q = [q1; q2];
%     A = [A1 A2];
%     C = blkdiag(C1,C2);
%     d = [d1; d2];
%     sol = quadprog(P,q,-C,-d,A,b);
%     y1 = sol(1:end/2); y2 = sol(end/2+1:end);
    
    U1 = y1((N+1)*nx+1:(N+1)*nx+1+nu-1,1); U2 = y2((N+1)*nx+1:(N+1)*nx+1+nu-1,1);
    x1 = Ad(1,1)*x1 + Bd(1,1)*U1 + Bd(1,2)*U2; % agent 1
    x2 = Ad(2,2)*x2 + Bd(2,1)*U1 + Bd(2,2)*U2; % agent 2

    implementedU1(:,i) = U1;
    implementedU2(:,i) = U2;
    xs1(:,i+1) = x1;
    xs2(:,i+1) = x2;
    
    
%     figure(1)
%     plot(times,xs1(:))
%     hold on
%     plot(times,xs2(:))
%     hold off
    i
    niters = [niters info.iters];
end


% Need to turn continuous-time problem into discrete-time problem next
%% View iterations requierd for each optimisation problem
figure(2)
plot(times(1:end-1),niters,'b','LineWidth',1.5);
hold on

%% Run simulation with warmstart

alpha0 = 0.1; gamma = 1; tau = 1e-4; show = false; visualise = false;
alpha0 = 0.45;
Ts = 0.5;
Tfinal = 20;
nsteps = round(Tfinal/Ts);
times = 0:Ts:Tfinal;

x1 = 80;
x2 = 20;
xs1 = zeros(nx,nsteps+1);
xs1(:,1) = x1;
xs2 = zeros(nx,nsteps+1);
xs2(:,1) = x2;

implementedU1 = zeros(nu,nsteps);
implementedU2 = zeros(nu,nsteps);
niters = [];
for i = 1:nsteps
    [P1,P2,q1,q2,A1,A2,b,C1,d1,C2,d2] = mpc2dualdecomp(Amat,Bmat,Q1,Q2,R1,R2,x1set,x2set,u1min,u1max,u2min,u2max,x1min,x1max,x2min,x2max,N,x1,x2);
    
    if (i > 1)
        [y1,y2,info] = dual_decomp(P1,P2,q1,q2,A1,A2,b,C1,d1,C2,d2,alpha0,gamma,tau,show,visualise,info);
    else
        [y1,y2,info] = dual_decomp(P1,P2,q1,q2,A1,A2,b,C1,d1,C2,d2,alpha0,gamma,tau,show,visualise);
    end
%     P = blkdiag(P1,P2);
%     q = [q1; q2];
%     A = [A1 A2];
%     C = blkdiag(C1,C2);
%     d = [d1; d2];
%     sol = quadprog(P,q,-C,-d,A,b);
%     y1 = sol(1:end/2); y2 = sol(end/2+1:end);
    
    U1 = y1((N+1)*nx+1:(N+1)*nx+1+nu-1,1); U2 = y2((N+1)*nx+1:(N+1)*nx+1+nu-1,1);
    x1 = Ad(1,1)*x1 + Bd(1,1)*U1 + Bd(1,2)*U2; % agent 1
    x2 = Ad(2,2)*x2 + Bd(2,1)*U1 + Bd(2,2)*U2; % agent 2

    implementedU1(:,i) = U1;
    implementedU2(:,i) = U2;
    xs1(:,i+1) = x1;
    xs2(:,i+1) = x2;
    
    
%     figure(1)
%     plot(times,xs1(:))
%     hold on
%     plot(times,xs2(:))
%     hold off
    i
    
    niters = [niters info.iters];

end
%% View iterations requierd for each optimisation problem
figure(2)
plot(times(1:end-1),niters,'r','LineWidth',1.5);
%% Figure labels
xlabel('Time (s)')
ylabel('Number of Iterations')
legend('No warmstart','Warmstart')