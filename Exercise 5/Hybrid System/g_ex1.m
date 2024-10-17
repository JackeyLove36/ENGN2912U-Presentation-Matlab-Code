function xplus = g_ex1(x)
% DJscription: jump map
% state
x1 = x(1);
x2 = x(2);
lambda = 0.8;
xplus = [-x1 ; -lambda*x2];
end