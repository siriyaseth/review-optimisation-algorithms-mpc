function [dotxi] = tilted_table_ode(t,xi,mu)
%TILTED_TABLE_ODE 
%   ode for tilted table and ball toy
g = -9.81;
theta = mu(1); phi = mu(2);
dotx = xi(2); doty = xi(4);
ddotx = 5.0/7.0*g*sin(theta);
ddoty = 5.0/7.0*g*sin(phi);

dotxi = [dotx; ddotx; doty; ddoty];
end

