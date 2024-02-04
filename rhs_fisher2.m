%Autores: Unai Gurbindo y Aitor Ayape
%Grado: Ciencia de Datos
%Modelización y Simulación de Sistemas Biológicos
% Cuestión 2 - Proyecto 4.1
%__________________________________________________________________________
function udot = rhs_fisher2(t,tipo,u)
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
% Coeficiente de difusión, D:
if tipo == 2
    D = 0;
else
    D = 1;
end
%__________________________________________________________________________
% Coeficiente de reacción, k = k(t):
if tipo==1
    k = 0;
elseif tipo==2
    k = 0.65;
elseif tipo==3
    k = 0.65;
elseif tipo==4
    k = 0.1625*(1-tanh(t-5))^2;
else
    k = -0.65;
end 
%__________________________________________________________________________
%Tipos
%1: difusión pura
%2: reacción pura
%3: reacción-difusión con generación ilimitada de biomasa
%4: reacción-difusión con generación limitada de biomasa
%5: reacción-difusión con destrucción de biomasa
%__________________________________________________________________________
% Sistema de EDOs:
udot = (1/h^2)*D*A*u + k*u.*(1-u);
end