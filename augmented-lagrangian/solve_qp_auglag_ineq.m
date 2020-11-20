function [xout,k] = solve_qp_auglag_ineq(P,q,A,b,C,d,mu0,gamma1,tau,alpha0,gamma2,show,visualise)
%SOLVE_QP_AUGLAG_BOUND
%   Solve QP problem with inequality constraints using augmented lagrangian

nI = size(C,1); %number of inequality constraints
nX = size(P,1); % number of decision variables in origianl problem
nE = size(A,1); % number od equality constraints in original problem
x0 = zeros(nX,1);
s0 = zeros(nI,1);

k = 1;
x = x0; s = s0;
y = [x;s];
mu = mu0;

lambdaE = zeros(nE,1);
lambdaI = zeros(nI,1);
lambda = [lambdaE; lambdaI];

primal = [];
dual = [];
while (true)   
    Ps = blkdiag(P, zeros(nI)); qs = [q; zeros(nI,1)];
    As = [A zeros(nE,nI);
          C -eye(nI)];
    bs = [b; d];
%     Cs = [zeros(nI,nX) eye(nI)];
%     ds = zeros(nI,1);

    P_La = Ps + mu*transpose(As)*As;
    q_La = qs - transpose(As)*lambda(:,k) - mu*transpose(As)*bs;
%     test_sol = quadprog(P_La,q_La,-Cs,-ds);
    
    l = [-inf(nX,1); zeros(nI,1)];
    u = inf(nX+nI,1);
    y_solved = solve_qp_grad_proj_box(P_La, q_La, l, u, alpha0, gamma2, tau, false, zeros(nX+nI,1));
%     y_solved = solve_qp_grad_proj_box(P_La, q_La, l, u, alpha0, gamma2, tau, true, zeros(nX+nI,1));
    y(:,k+1) = y_solved(:,end);
    
    x(:,k+1) = y(1:nX,k+1);
    s(:,k+1) = y(nX+1:end,k+1);

    lambda(:,k+1) = lambda(:,k) - mu*(As*y(:,k+1) - bs);
    lambdaE(:,k+1) = lambda(1:nE,k+1);
    lambdaI(:,k+1) = lambda(nE+1:end,k+1);
    
    if (show)
        primal = [primal primal_fun(x(:,k+1),P,q)];
        dual = [dual dual_fun(x(:,k+1), lambdaE(:,k+1), lambdaI(:,k+1), P, q, A, b, C, d)];
        fprintf('k=%d: f(x)=%.16f, L(x,lambdaE,lambdaI)=%.16f\n', k, primal(k), dual(k));
    end

    if (norm(lambda(:,k+1)-lambda(:,k),2) < tau) % q(mu,lambda1,lambda2) concave function
        break
    end
    k = k + 1;
    mu = mu*gamma1;
    if (visualise)
        figure(1)
        plot(primal);
        hold on
        plot(dual);
        hold off
        pause(0.001)
    end

end
xout = x(:,end);
end

% Helper functions

function p = primal_fun(x,P,q)
    p = 0.5*transpose(x)*P*x + transpose(q)*x;
end

function q = dual_fun(x,lambdaE,lambdaI,P,q,A,b,C,d)
    q = primal_fun(x,P,q) - transpose(lambdaE)*(A*x - b) - transpose(lambdaI)*(C*x-d);
end

