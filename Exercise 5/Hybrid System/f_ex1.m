function xdot = f_ex1(x)
% Description: Flow map   
% state
x2 = x(2);

% differential equations
gamma = -9.81;
xdot = [x2 ; gamma];
end