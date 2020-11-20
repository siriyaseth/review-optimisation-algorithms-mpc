% Function for solving QP problem with box constraints using gradient
% projection and diminishing step length but no armijo rule
function [x,k] = solve_qp_grad_proj_box(Q, q, l, u, alpha0, gamma, tau, show, x0)
if nargin == 8
    x0 = zeros(size(Q,1),1);
end
x = x0;
d = [];
m = [];
if (show)
    fprintf('k = %d: f(x) = %.16f\n', 0, objective(Q,q,x(:,1)));
end
k = 1;
alpha = alpha0;
while (true)
    x(:,k+1) = project_box(x(:,k) - alpha*gradfunc(Q,q,x(:,k)),l,u);
    if (show)
        fprintf('k = %d: f(x) = %.16f, ||y(k+1)-y(k)|| = %.16f\n', k, objective(Q,q,x(:,k+1)), norm(x(:,k+1) - x(:,k),2));
    end
    if ( norm(x(:,k+1) - x(:,k),2) < tau)
        break;
    end
    if (k > 10000)
        % fail to converge
        break;
    end
    k = k + 1;
    alpha = alpha*gamma;
end
end

function fx = objective(Q,q,x)
fx = 0.5*x'*Q*x + q'*x;
end

function gradfx = gradfunc(Q,q,x)
gradfx = Q*x + q;
end

% Projection for general box constraint
function p = project_box(x,l,u)
%     p =  min([max([x,l]),u]);
    p = min(max(x,l),u);
end


