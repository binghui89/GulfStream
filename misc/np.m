function np()
fun = @DIST; % Distance function
x0 = [0; 0; 1];
A = [];
b = [];
Aeq = [];
beq = [];
lb = [];
ub = [];
nonlcon = @mycon;
x = fmincon(fun,x0,A,b,Aeq,beq,lb,ub,nonlcon);
fprintf('a = %f, b = %f, c = %f \n', x(1), x(2), x(3));
end

% x = [a, b, c]'
function F = DIST(x)
cap = [4, 40, 200, 400]; % MW
cost = [10783335, 22271891, 42697218, 48954903]; % $
F = ( cost(1) - x(1) - x(2)*log(cap(1)+x(3)) )^2 +...
    ( cost(2) - x(1) - x(2)*log(cap(2)+x(3)) )^2 +...
    ( cost(3) - x(1) - x(2)*log(cap(3)+x(3)) )^2 +...
    ( cost(4) - x(1) - x(2)*log(cap(4)+x(3)) )^2;
end

function [c,ceq, gradc, gradceq] = mycon(x)
c = [];
ceq = x(1) + x(2)*log(x(3));
gradc = [];
gradceq = [1; log(x(3)); x(2)/x(3)];
end