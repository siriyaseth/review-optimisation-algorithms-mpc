function [yout1,yout2,info] = dual_decomp(P1,P2,q1,q2,A1,A2,b,C1,d1,C2,d2,alpha0,gamma,tau,show,visualise,warmstart)
%DUAL_DECOMP
%   Code for running dual decomosition based on problem setup
    
if (nargin < 17)
    y10 = zeros(size(P1,1),1); y20 = zeros(size(P2,1),1);
    mu0 = zeros(size(A1,1),1);
    lambda10 = zeros(size(C1,1),1); lambda20 = zeros(size(C2,1),1); 
else
    y10 = warmstart.y1;
    y20 = warmstart.y2;
    mu0 = warmstart.mu;
    lambda10 = warmstart.lambda1;
    lambda20 = warmstart.lambda2;
end



k = 1;
y1 = y10; y2 = y20;
alpha = alpha0;

mu = mu0; lambda1 = lambda10; lambda2 = lambda20;
primal = []; dual = [];

while (true)

    lagrange_q1 = (q1 - A1'*mu(:,k) - C1'*lambda1(:,k));
    lagrange_q2 = (q2 - A2'*mu(:,k) - C2'*lambda2(:,k));
    y1(:,k+1) = P1\(-lagrange_q1);
    y2(:,k+1) = P2\(-lagrange_q2);
    mu(:,k+1) = mu(:,k) - alpha*(A1*y1(:,k) + A2*y2(:,k) - b);
    lambda1(:,k+1) = proj_pos(lambda1(:,k) - alpha*(C1*y1(:,k) - d1));
    lambda2(:,k+1) = proj_pos(lambda2(:,k) - alpha*(C2*y2(:,k) - d2));
    primal = [primal primal_fun(y1(:,k+1),y2(:,k+1),P1,P2,q1,q2)];
    dual = [dual dual_fun(y1(:,k+1),y2(:,k+1),mu(:,k+1),lambda1(:,k+1),lambda2(:,k+1),P1,P2,q1,q2,A1,A2,b,C1,d1,C2,d2)];
    if (show)
        fprintf('k=%d: f(y1,y2)=%.16f, q(mu,lambda1,lambda2)=%.16f\n', k, primal(k), dual(k));
    end

    if (norm([mu(:,k) - mu(:,k+1); lambda1(:,k) - lambda1(:,k+1); lambda2(:,k) - lambda2(:,k+1)],2) < tau) % q(mu,lambda1,lambda2) concave function
        break
    end
    k = k + 1;
    alpha = alpha*gamma;
    
    if (visualise)
        figure(1)
        plot(primal);
        hold on
        plot(dual);
        hold off
        pause(0.001)
    end
end

yout1 = y1(:,end); yout2 = y2(:,end);

info = struct;
info.y1 = yout1;
info.y2 = yout2;
info.mu = mu(:,end);
info.lambda1 = lambda1(:,end);
info.lambda2 = lambda2(:,end);
info.iters = k;

end

% Project onto nonnegative orthant
function p = proj_pos(x)
n = size(x,1);
p = max(zeros(n,1),x);
end

function p = primal_fun(y1,y2,P1,P2,q1,q2)
    p = 0.5*(transpose(y1)*P1*y1 + transpose(y2)*P2*y2) + transpose(q1)*y1 + transpose(q2)*y2;
end

function q = dual_fun(y1,y2,mu,lambda1,lambda2,P1,P2,q1,q2,A1,A2,b,C1,d1,C2,d2)
    q = primal_fun(y1,y2,P1,P2,q1,q2) - transpose(mu)*(A1*y1 + A2*y2 - b) - transpose(lambda1)*(C1*y1 - d1) - transpose(lambda2)*(C2*y2 - d2);
end
