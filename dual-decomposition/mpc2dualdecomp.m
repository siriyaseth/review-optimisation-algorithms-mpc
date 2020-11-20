function [P1,P2,q1,q2,A1,A2,b,C1,d1,C2,d2] = mpc2dualdecomp(A,B,Q1,Q2,R1,R2,x1set,x2set,u1min,u1max,u2min,u2max,x1min,x1max,x2min,x2max,N,x10,x20)
%MPC2DUALDECOMP
%   Converts 2-agent distributed MPC problem into dual decomposition
%   optimization problem

nx = size(B{1,1},1); nu = size(B{1,1},2);

Q1bar = kron(eye(N+1), Q1); R1bar = kron(eye(N), R1);
Q2bar = kron(eye(N+1), Q2); R2bar = kron(eye(N), R2);

P1 = blkdiag(Q1bar, R1bar); P2 = blkdiag(Q2bar, R2bar);

y1set = [kron(ones(N+1,1),x1set); zeros(N*nu,1)];
y2set = [kron(ones(N+1,1),x2set); zeros(N*nu,1)];

q1 = -P1*y1set; q2 = -P2*y2set;

D11 = [zeros(nx, (N+1)*nx + N*nu); 
       kron(eye(N), A{1,1}) zeros(N, nx) kron(eye(N), B{1,1})];

D12 = [zeros(nx, (N+1)*nx + N*nu); 
       kron(eye(N), A{1,2}) zeros(N, nx) kron(eye(N), B{1,2})];

D21 = [zeros(nx, (N+1)*nx + N*nu); 
       kron(eye(N), A{2,1}) zeros(N, nx) kron(eye(N), B{2,1})];

D22 = [zeros(nx, (N+1)*nx + N*nu); 
       kron(eye(N), A{2,2}) zeros(N, nx) kron(eye(N), B{2,2})];
   
ny = size(y1set,1);
nxs = (N+1)*nx;
nus = N*nu;
A1 = [[eye(nxs) zeros(nxs,nus)] - D11; -D21]; A2 = [-D12; [eye(nxs) zeros(nxs,nus)] - D22];

xs10 = [x10; zeros(N*nx, 1)];
xs20 = [x20; zeros(N*nx, 1)];
b = [xs10; xs20];

y1min = [kron(ones(N+1,1),x1min); kron(ones(N,1),u1min)];
y1max = [kron(ones(N+1,1),x1max); kron(ones(N,1),u1max)];
y2min = [kron(ones(N+1,1),x2min); kron(ones(N,1),u2min)];
y2max = [kron(ones(N+1,1),x2max); kron(ones(N,1),u2max)];

C1 = [eye(ny); -eye(ny)];
d1 = -[-y1min; y1max];
C2 = [eye(ny); -eye(ny)];
d2 = -[-y2min; y2max];

% Modificationt so that first state has no inequality constraints applied
C1 = [blkdiag(0,eye(ny-1)); -blkdiag(0,eye(ny-1))];
C2 = [blkdiag(0,eye(ny-1)); -blkdiag(0,eye(ny-1))];


end