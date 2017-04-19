% Chebyshev class test script

clear all; clc

% initialization
A = Chebyshev(3,'x');   % T_3(x)
B = Chebyshev(2,'x');   % T_2(x)
C = Chebyshev(2,'y');   % T_2(y)
D = Chebyshev(3,'y');   % T_3(y)
F = Chebyshev(3,'z');   % T_3(z)
G = Chebyshev(4,'x');   % T_4(x)

% copy a polynomial
J = Chebyshev(A)    % T_3(x)

% polynomial using list of coefficients with descending degree, similar to poly()
K = Chebyshev([1 2 3 4],'x')    % T_3(x) + 2*T_2(x) + 3*T_1(x) + 4

% test basic multiplication
X = A*B         % T_3(x)*T_2(x)
Y = B*C         % T_2(x)*T_2(y)

% test distributivity and addition
Z = C*(A+B)     % T_2(y)*(T_3(x) + T_2(x)) = T_2(y)*T_3(x) + T_2(y)*T_2(x)

% coefficients with distributivity
Q = 2*C*(A+B)   % 2*T_2(y)*(T_3(x) + T_2(x))) = 2*T_2(y)*T_3(x) + 2*T_2(y)*T_2(x)

% coefficient multiplication
P = 3*2*A       % 3*2*T_3(x) = 6*T_3(x)

% test sort()
R = A*C*B*D     % T_3(x)*T_2(y)*T_2(x)*T_3(y)

% groups product in terms of argument variable
sort(R)         % = T_3(x)*T_2(x)*T_2(y)*T_3(y)

% multi-variable
T = A*C*B*D*F   % T_3(x)*T_2(y)*T_2(x)*T_3(y)*T_3(z)

% differentiation
diff(A,'x')     % 6*T_1(x) + 6*T_3(x) + 6*T_5(x)

% test differentiation: product rule
W = diff(X,'x')
% diff(X) = diff(T_3(x))*T_2(x) + T_3(x)*diff(T_2(x))
%         = 6*T_1(x)*T_2(x) + 6*T_3(x)*T_2(x) + 6*T_5(x)*T_2(x) + 4*T_1(x)*T_3(x) + 4*T_3(x)*T_3(x)

% test simplification
simplify(W)     % 3*T_3(x) + 3*T_1(x) + 3*T_5(x) + 3*T_1(x) + 3*T_7(x) + 3*T_3(x) + 2*T_4(x) + 2*T_2(x) + 2*T_6(x) + 2

% test differentiation: sum
S = A + B       % T_3(x) + T_2(x)
H = diff(S,'x') % 6*T_1(x) + 6*T_3(x) + 6*T_5(x) + 4*T_1(x) + 4*T_3(x)
