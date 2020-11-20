function [P,q] = mpc2qp_compact(x0,A,B,Q,R,N)
%MPC2QP_INPUT Optains objective function in QP form from MPC in compact
%form
% 
% x0 = [3;2];
% A = ones(2,2); B = ones(2,3);
% N = 4;
% Q = eye(2);
% R = 4

Qbar = kron(eye(N), Q);
Rbar = kron(eye(N), R);

nx = size(B,1);
nu = size(B,2);
G = zeros(N*nx, N*nu);
for j = 1:N
    cntr = 0;
    for i = j:N
        G((i-1)*nx+1:i*nx, (j-1)*nu+1:j*nu) = A^cntr*B;
        cntr = cntr + 1;
    end
end

c = zeros(N*nx, 1);
for i = 1:N
    c((i-1)*nx+1:i*nx, 1) = A^i * x0;
end

P = G'*Qbar*G + Rbar;
q = G'*Qbar*c;
end

