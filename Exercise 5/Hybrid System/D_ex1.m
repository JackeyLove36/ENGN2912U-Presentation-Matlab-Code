function inside = D_ex1(x) 
% Description: Jump set
x1 = x(1);
x2 = x(2);
if (x1 <= 0 && x2 <= 0)
    inside = 1;
else
    inside = 0;
end
end