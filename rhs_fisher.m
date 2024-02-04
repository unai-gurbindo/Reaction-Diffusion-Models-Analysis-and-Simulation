%Autores: Unai Gurbindo y Aitor Ayape
%Grado: Ciencia de Datos
%Modelizaci�n y Simulaci�n de Sistemas Biol�gicos
% Cuesti�n 2 - Proyecto 4.1
%__________________________________________________________________________
function udot = rhs_fisher(t,u)
%__________________________________________________________________________
%DATOS
global h N
%__________________________________________________________________________
% Matriz de coeficientes, A:
e = ones(N+2,1);
C = spdiags([e -2*e e],-1:1,N+2,N+2);
%Completamos la matriz C
C(1,2)=2;
C(N+2,N+1)=2;
I = speye(N+2);
A = kron(I,C)+kron(C,I);
%__________________________________________________________________________
% Coeficiente de difusi�n, D:
%D = 0;
D = 1;
%__________________________________________________________________________
% Coeficiente de reacci�n, k = k(t):
%k = 0;
%k = 0.65;
%k = 0.1625*(1-tanh(t-5))^2;
k = -0.65;
%__________________________________________________________________________
% Sistema de EDOs:
udot = (1/h^2)*D*A*u + k*u.*(1-u);

end