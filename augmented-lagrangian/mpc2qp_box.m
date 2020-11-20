function [P,q,A,b,C,d] = mpc2qp_box(Phi,Gamma,Q,R,xset,umin,umax,xmin,xmax,N,x0)
%MPC2DUALDECOMP
%   Converts 2-agent distributed MPC problem into dual decomposition
%   optimization problem


nx = size(Gamma,1);
nu = size(Gamma,2);

% Setup matrices for box constraints on states/inputs
F = [eye(nx);
     -eye(nx);
     zeros(2*nu,nx)];
G = [zeros(2*nx,nu);
     eye(nu);
     -eye(nu)];
beta = [xmin; -xmax; umin; -umax];
H = [eye(nx);
     -eye(nx)];
kappa = [xmin; -xmax];


P = blkdiag(kron(eye(N+1),Q),kron(eye(N),R));
Xbar = kron(ones(N+1,1), xset);
ybar = [Xbar; zeros(N*nu,1)];

q = -P*ybar;

D = [zeros(nx,N*nx) zeros(nx,nx) zeros(nx, N*nu);
     kron(eye(N),Phi) zeros(N*nx,nx) kron(eye(N),Gamma)];
A = [eye((N+1)*nx) zeros((N+1)*nx,N*nu)] - D;

b = [x0; zeros(N*nx,1)];

C = [kron(eye(N),F) zeros(N*size(F,1),size(H,2)) kron(eye(N),G);
     zeros(size(H,1),N*size(F,2)) H zeros(size(H,1),N*size(G,2))];

d = [kron(ones(N,1),beta); kappa];

% Ensure inequality constraints are not enforced at first state
C = C(nx+1:end,1:end);
d = d(nx+1:end,1);

end