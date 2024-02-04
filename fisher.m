%Autores: Unai Gurbindo y Aitor Ayape
%Grado: Ciencia de Datos
%Modelización y Simulación de Sistemas Biológicos
% Cuestión 2 - Proyecto 4.1
%__________________________________________________________________________
% Modelo de dispersión genética (ecuación de Fisher):
%   u_t = D*(u_xx+u_yy)+k*u*(1-u).
% Discretización mediante un esquema de diferencias finitas de
% 5 puntos en espacio combinado con diferentes integradores temporales
% (método de líneas).
%__________________________________________________________________________
%IMPLEMENTACIÓN:
%--------------------------------------------------------------------------
clear all
%__________________________________________________________________________
%DATOS DEL PROBLEMA
global h N

v = VideoWriter('fisher.avi');
v.FrameRate = 4;
open(v)

a = -2 ;
b = 2 ;

N = 47;             % Número de nodos interiores en las direcciones x e y
h = (b-a)/(N+1) ;             % Tamaño de malla
 
x = linspace(a,b,N+2);
y = linspace(a,b,N+2);

t0 = 0 ;
tf = 10;
%__________________________________________________________________________
% Condición inicial:
[X,Y] = meshgrid(x,y);

r = 0.35 ;
x0 = 0.5 ; y0 = -0.5 ;
x1 = -0.5 ; y1 = 0.5 ;
Z1 = (1/(2*pi*r^2)) * exp(-((X-x0).^2+(Y-y0).^2)/(2*r^2)) ;
Z2 = (1/(2*pi*r^2)) * exp(-((X-x1).^2+(Y-y1).^2)/(2*r^2)) ;
u0 = (Z1+Z2)';

u0 = reshape(u0,(N+2)*(N+2),1);
%__________________________________________________________________________
% Integración en tiempo con ode45/ode15s:
options = odeset('RelTol',1.e-6,'AbsTol',1.e-6,'Stats','on','Refine',1);
%__________________________________________________________________________
% Rutina ode45/ode15s:
fprintf('\n')
fprintf('Estadísticas\n');
fprintf('------------\n');

%[t,u] = ode45(@rhs_fisher,[t0,tf],u0,options); %ode45
[t,u] = ode15s(@rhs_fisher,[t0,tf],u0,options); %ode15s

% Representación gráfica del tamaño de paso en tiempo: 
% figure(1)
% semilogy(t(2:end),diff(t),'o','Markersize',5);
% xlabel('$t$','Interpreter','latex')
% ylabel('$\Delta t$','Interpreter','latex')
% title('Tamaño de paso en tiempo')

% % Integración en tiempo con Euler explícito:
% M = 1450;
% u = Euler('rhs_fisher',[t0,tf],u0,M);
% t = u(:,1);
% u(:,1) = [];

for i = 1:size(t)
    sol = reshape(u(i,:),N+2,N+2);
    figure(2)
    surf(x,y,sol');
    %contourf(x,y,sol');
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$u(x,y)$','Interpreter','latex');
    title(['Solución numérica, N = ',num2str(N),', t = ',num2str(t(i),'%6.4f')]);
    %title(['Solución numérica, N = ',num2str(N),', M = ',num2str(M),', t = ',num2str(t(i),'%6.4f')]);
    colorbar;
    colormap('jet');
    set(gca,'CLim',[0 1.3]);
    axis([a b a b 0 1.3]);
    frame=getframe(gcf);
    writeVideo(v,frame);
end

close(v);