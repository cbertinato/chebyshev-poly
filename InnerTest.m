% Testing inner()
clear A B C D E X Y; clc

A = Chebyshev(2,'x');
B = Chebyshev(3,'x');
C = Chebyshev(2,'y');
D = Chebyshev(3,'y');
E = Chebyshev(0,'x');

X = A*B;
fprintf('inner(T_2(x)*T_3(x)), Expected: 0\n')
inner(X)

Y = A*A;
fprintf('inner(T_2(x)*T_2(x)), Expected: pi/2 = 1.57079...\n')
inner(Y)

Y = 2*A*A;
fprintf('inner(2*T_2(x)*T_2(x)), Expected: pi\n')
inner(Y)

Y = E*E;
fprintf('inner(T_0(x)*T_0(x)), Expected: pi\n')
inner(Y)

Y = A*A*C;
fprintf('inner(T_2(x)*T_2(x)*T_2(y)), Expected: 0\n')
inner(Y)

Y = 3*A*B*C*C;
fprintf('inner(3*T_2(x)*T_3(x)*T_2(y)*T_2(y)), Expected: 0\n')
inner(Y)

Y = 3*A*A*C*C;
fprintf('inner(3*T_2(x)*T_2(x)*T_2(y)*T_2(y)), Expected: 3*pi^2/4 = 7.4022...\n')
inner(Y)

Y = (A+B)*(C+D);
fprintf('inner(T_2(x)*T_2(y) + T_2(x)*T_3(y) + T_3(x)*T_2(y) + T_3(x)*T_3(y)), Expected: 0\n')
inner(Y)

fprintf('inner(T_2(x) + T_3(x),T_2(x) + T_3(x)), Expected: pi\n')
inner(A+B,A+B)

Y = (A+B)*(A+B);
fprintf('inner(T_2(x)*T_2(x) + T_2(x)*T_3(x) + T_3(x)*T_2(x) + T_3(x)*T_3(x)), Expected: pi\n')
inner(Y)