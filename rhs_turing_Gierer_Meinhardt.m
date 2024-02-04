%Autores: Unai Gurbindo y Aitor Ayape
%Grado: Ciencia de Datos
%Modelización y Simulación de Sistemas Biológicos
% Cuestiones 1 y 2 - Proyecto 4.2
%__________________________________________________________________________
function wdot = rhs_turing(t,w)
%__________________________________________________________________________
%DATOS
global a b N hx hy
%__________________________________________________________________________
% Parámetros del modelo:
Du = 0.05;

Dv = 1;
gamma = 100;

u = w(1:(N+2)^2);
v = w((N+2)^2+1:2*(N+2)^2);
%__________________________________________________________________________
% Matriz de coeficientes, M:
%--------------------------------------------------------------------------
% Modelo de Gierer-Meinhardt:
f = gamma*(a-b*u+u.^2./v);
g = gamma*(u.^2-v);
%--------------------------------------------------------------------------
% Matriz de coeficientes, A:
e = ones(N+2,1);
C = spdiags([e -2*e e],-1:1,N+2,N+2);
C(1,2) = 2;
C(end,end-1) = 2;

I1 = speye(N+2)./hx^2;
I2 = speye(N+2)./hy^2;
A = kron(I1,C)+kron(C,I2);
%--------------------------------------------------------------------------
M1=Du.*A;
M2=Dv.*A;
M= blkdiag(M1,M2);
%__________________________________________________________________________
% Sistema de EDOs:
wdot = M*w+[f;g];

end