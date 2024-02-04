%Autores: Unai Gurbindo y Aitor Ayape
%Grado: Ciencia de Datos
%Modelización y Simulación de Sistemas Biológicos
% Cuestiones 1 y 2 - Proyecto 4.2
%__________________________________________________________________________
% Modelo de formación de patrones (sistema de Turing):
%   u_t = Du*(u_xx+u_yy)+f(u,v),
%   v_t = Dv*(v_xx+v_yy)+g(u,v).
% Discretización mediante un esquema de diferencias finitas de
% 5 puntos en espacio combinado con diferentes integradores temporales
% (método de líneas).
%__________________________________________________________________________
%__________________________________________________________________________
%IMPLEMENTACIÓN:
%--------------------------------------------------------------------------
clear all;
close all;
diary on;
mkdir("Evaluaciones-4.2-C1-C2");
%Gierer-Meinhardt
mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt");
mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Cuadrado");
mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Cuadrado/Rutina-ode45");
mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Cuadrado/Rutina-ode15s");
mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Cuadrado/Rutina-ode45/Componente-U");
mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Cuadrado/Rutina-ode15s/Componente-U");
mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Cuadrado/Rutina-ode45/Componente-V");
mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Cuadrado/Rutina-ode15s/Componente-V");
mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Rectangulo");
% mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Rectangulo/Rutina-ode45");
mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Rectangulo/Rutina-ode15s");
% mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Rectangulo/Rutina-ode45/Componente-U");
mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Rectangulo/Rutina-ode15s/Componente-U");
% mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Rectangulo/Rutina-ode45/Componente-V");
mkdir("./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Rectangulo/Rutina-ode15s/Componente-V");
%Schnakenberg
mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg");
mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Cuadrado");
mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Cuadrado/Rutina-ode45");
mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Cuadrado/Rutina-ode15s");
mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Cuadrado/Rutina-ode45/Componente-U");
mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Cuadrado/Rutina-ode15s/Componente-U");
mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Cuadrado/Rutina-ode45/Componente-V");
mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Cuadrado/Rutina-ode15s/Componente-V");
mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Rectangulo");
% mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Rectangulo/Rutina-ode45");
mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Rectangulo/Rutina-ode15s");
% mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Rectangulo/Rutina-ode45/Componente-U");
mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Rectangulo/Rutina-ode15s/Componente-U");
% mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Rectangulo/Rutina-ode45/Componente-V");
mkdir("./Evaluaciones-4.2-C1-C2/Schnakenberg/Rectangulo/Rutina-ode15s/Componente-V");
%__________________________________________________________________________
%DATOS Cuadrado
global a b N hx hy
% Parámetros del modelo:
a = 0.1305 ;
b = 0.7695 ;
% Parámetros de la discretización:
ax = 0;
bx = 1;
ay = 0;
by = 1;
N = 47;% Número de nodos interiores en las direcciones x e y
hx = (bx-ax)/(N+1);      % Tamaño de malla en la dirección x
hy = (by-ay)/(N+1);      % Tamaño de malla en la dirección y 
x = linspace(ax,bx,N+2);
y = linspace(ay,by,N+2);
t0 = 0;
[X,Y] = meshgrid(x,y);% Condición inicial
options = odeset('RelTol',1.e-6,'AbsTol',1.e-6,'Stats','on','Refine',1);% Integración en tiempo
%__________________________________________________________________________
%TIEMPOS DE EVALUACIÓN
tiempos=[0.5 1 2 10];
%__________________________________________________________________________
%Cuadrado
for temp = 1:4
    tf=tiempos(temp);
    %----------------------------------------------------------------------
    %Gierer-Meinhardt
    u_star = (a+1)/b;
    v_star = u_star^2;
    u0 = u_star + 10^-3 * exp(-100*((X-(bx-ax)/3).^2+(Y-(by-ay)/3).^2));
    v0 = v_star * ones(N+2);    
    u0 = reshape(u0,(N+2)*(N+2),1);
    v0 = reshape(v0,(N+2)*(N+2),1);
    w0 = [u0;v0];
    %......................................................................
    %Rutina ode45
    fprintf('\n')
    fprintf(sprintf('Estadísticas, ode45, Gierer-Meinhardt, tiempo= %d, Cuadrado\n',tf));
    fprintf('------------\n');
    [t,w] = ode45(@rhs_turing_Gierer_Meinhardt,[t0,tf],w0,options);
    usol = reshape(w(length(t),1:(N+2)^2),N+2,N+2);
    vsol = reshape(w(length(t),(N+2)^2+1:2*(N+2)^2),N+2,N+2);
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %Cuadrado
    %Componente U
    figure;
    %     surf(x,y,usol');
    contourf(x,y,usol');
    colorbar;    
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$u(x,y)$','Interpreter','latex');
    title(['Solución numérica, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
    set(gca,'CLim',[-1.5 5]);
    axis([ax bx ay by -1.5 5]);
    fileName=sprintf('./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Cuadrado/Rutina-ode45/Componente-U/Tiempo-%d-turing-contourf.png',tf);
    saveas(gcf,fileName);
    close all;
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %Componente V      
    figure;
    %     surf(x,y,vsol');
    contourf(x,y,vsol');
    colorbar;    
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$v(x,y)$','Interpreter','latex');
    title(['Solución numérica, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
    set(gca,'CLim',[-1.5 5]);
    axis([ax bx ay by -1.5 5]);
    fileName=sprintf('./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Cuadrado/Rutina-ode45/Componente-V/Tiempo-%d-turing-contourf.png',tf);
    saveas(gcf,fileName);
    %......................................................................
    %Rutina ode15s
    fprintf('\n')
    fprintf(sprintf('Estadísticas, ode15s, Gierer-Meinhardt, tiempo= %d, Cuadrado\n',tf));
    fprintf('------------\n');
    [t,w] = ode15s(@rhs_turing_Gierer_Meinhardt,[t0,tf],w0,options);
    usol = reshape(w(length(t),1:(N+2)^2),N+2,N+2);
    vsol = reshape(w(length(t),(N+2)^2+1:2*(N+2)^2),N+2,N+2);
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %Cuadrado
    %Componente U
    figure;
    %     surf(x,y,usol');
    contourf(x,y,usol');
    colorbar;    
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$u(x,y)$','Interpreter','latex');
    title(['Solución numérica, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
    set(gca,'CLim',[-1.5 5]);
    axis([ax bx ay by -1.5 5]);
    fileName=sprintf('./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Cuadrado/Rutina-ode15s/Componente-U/Tiempo-%d-turing-contourf.png',tf);
    saveas(gcf,fileName);
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %Componente V      
    figure;
    %     surf(x,y,vsol');
    contourf(x,y,vsol');
    colorbar;    
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$v(x,y)$','Interpreter','latex');
    title(['Solución numérica, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
    set(gca,'CLim',[-1.5 5]);
    axis([ax bx ay by -1.5 5]);
    fileName=sprintf('./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Cuadrado/Rutina-ode15s/Componente-V/Tiempo-%d-turing-contourf.png',tf);
    saveas(gcf,fileName);
    %......................................................................
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %Schnakenberg
    u_star = (a+b);
    v_star = b/(u_star^2);
    u0 = u_star + 10^-3 * exp(-100*((X-(bx-ax)/3).^2+(Y-(by-ay)/3).^2));
    v0 = v_star * ones(N+2);    
    u0 = reshape(u0,(N+2)*(N+2),1);
    v0 = reshape(v0,(N+2)*(N+2),1);
    w0 = [u0;v0];
    %......................................................................
    %Rutina ode45
    fprintf('\n')
    fprintf(sprintf('Estadísticas, ode45, Schnakenberg, tiempo= %d, Cuadrado\n',tf));
    fprintf('------------\n');
    [t,w] = ode45(@rhs_turing_Schnakenberg,[t0,tf],w0,options);
    usol = reshape(w(length(t),1:(N+2)^2),N+2,N+2);
    vsol = reshape(w(length(t),(N+2)^2+1:2*(N+2)^2),N+2,N+2);
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %Cuadrado
    %Componente U
    figure;
    %     surf(x,y,usol');
    contourf(x,y,usol');
    colorbar;    
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$u(x,y)$','Interpreter','latex');
    title(['Solución numérica, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
    set(gca,'CLim',[-1.5 4]);
    axis([ax bx ay by -1.5 4]);
    fileName=sprintf('./Evaluaciones-4.2-C1-C2/Schnakenberg/Cuadrado/Rutina-ode45/Componente-U/Tiempo-%d-turing-contourf.png',tf);
    saveas(gcf,fileName);
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %Componente V      
    figure;
    %     surf(x,y,vsol');
    contourf(x,y,vsol');
    colorbar;    
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$v(x,y)$','Interpreter','latex');
    title(['Solución numérica, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
    set(gca,'CLim',[-1.5 4]);
    axis([ax bx ay by -1.5 4]);
    fileName=sprintf('./Evaluaciones-4.2-C1-C2/Schnakenberg/Cuadrado/Rutina-ode45/Componente-V/Tiempo-%d-turing.png',tf);
    saveas(gcf,fileName);
    %......................................................................
    %......................................................................
    %Rutina ode15s
    fprintf('\n')
    fprintf(sprintf('Estadísticas, ode15s, Schnakenberg, tiempo= %d, Cuadrado\n',tf));
    fprintf('------------\n');
    [t,w] = ode15s(@rhs_turing_Schnakenberg,[t0,tf],w0,options);
    usol = reshape(w(length(t),1:(N+2)^2),N+2,N+2);
    vsol = reshape(w(length(t),(N+2)^2+1:2*(N+2)^2),N+2,N+2);
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %Cuadrado
    %Componente U
    figure;
    %     surf(x,y,usol');
    contourf(x,y,usol');
    colorbar;    
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$u(x,y)$','Interpreter','latex');
    title(['Solución numérica, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
    set(gca,'CLim',[-1.5 4]);
    axis([ax bx ay by -1.5 4]);
    fileName=sprintf('./Evaluaciones-4.2-C1-C2/Schnakenberg/Cuadrado/Rutina-ode15s/Componente-U/Tiempo-%d-turing-contourf.png',tf);
    saveas(gcf,fileName);
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %Componente V      
    figure;
%     surf(x,y,vsol');
    contourf(x,y,vsol');
    colorbar;    
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$v(x,y)$','Interpreter','latex');
    title(['Solución numérica, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
    set(gca,'CLim',[-1.5 4]);
    axis([ax bx ay by -1.5 4]);
    fileName=sprintf('./Evaluaciones-4.2-C1-C2/Schnakenberg/Cuadrado/Rutina-ode15s/Componente-V/Tiempo-%d-turing-contourf.png',tf);
    saveas(gcf,fileName);
    %......................................................................
end
%__________________________________________________________________________
%DATOS Rectangulo
% Parámetros del modelo:
a = 0.1305 ;
b = 0.7695 ;
% Parámetros de la discretización:
ax = 0;
bx = 0.05;
ay = 0;
by = 1;
N = 47;% Número de nodos interiores en las direcciones x e y
hx = (bx-ax)/(N+1);      % Tamaño de malla en la dirección x
hy = (by-ay)/(N+1);      % Tamaño de malla en la dirección y 
x = linspace(ax,bx,N+2);
y = linspace(ay,by,N+2);
t0 = 0;
[X,Y] = meshgrid(x,y);% Condición inicial
%__________________________________________________________________________
%Rectangulo
for temp = 1:4
    tf=tiempos(temp);
    %----------------------------------------------------------------------
    %Gierer-Meinhardt
    u_star = (a+1)/b;
    v_star = u_star^2;
    u0 = u_star + 10^-3 * exp(-100*((X-(bx-ax)/3).^2+(Y-(by-ay)/3).^2));
    v0 = v_star * ones(N+2);    
    u0 = reshape(u0,(N+2)*(N+2),1);
    v0 = reshape(v0,(N+2)*(N+2),1);
    w0 = [u0;v0];
    %......................................................................
    %Rutina ode45
%     fprintf('\n')
%     fprintf(sprintf('Estadísticas, ode45, Gierer-Meinhardt, tiempo= %d, Rectangulo\n',tf));
%     fprintf('------------\n');
%     [t,w] = ode45(@rhs_turing_Gierer_Meinhardt,[t0,tf],w0,options);
%     usol = reshape(w(length(t),1:(N+2)^2),N+2,N+2);
%     vsol = reshape(w(length(t),(N+2)^2+1:2*(N+2)^2),N+2,N+2);
%     %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%     %Rectangulo
%     %Componente U
%     figure;
%     contourf(x,y,usol');
%     colorbar;    
%     daspect([1 3 1]);
%     xlabel('$x$','Interpreter','latex');
%     ylabel('$y$','Interpreter','latex');
%     zlabel('$u(x,y)$','Interpreter','latex');
%     title(['Solución numérica u, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
%     set(gca,'CLim',[0 4]);
%     axis([ax bx ay by 0 4]);
%     fileName=sprintf('./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Rectangulo/Rutina-ode45/Componente-U/Tiempo-%d-turing.png',tf);
%     saveas(gcf,fileName);
%     %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%     %Componente V      
%     figure;
%     contourf(x,y,vsol');
%     colorbar;    
%     daspect([1 3 1]);
%     xlabel('$x$','Interpreter','latex');
%     ylabel('$y$','Interpreter','latex');
%     zlabel('$v(x,y)$','Interpreter','latex');
%     title(['Solución numérica v, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
%     set(gca,'CLim',[0 4]);
%     axis([ax bx ay by 0 4]);
%     fileName=sprintf('./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Rectangulo/Rutina-ode45/Componente-V/Tiempo-%d-turing.png',tf);
%     saveas(gcf,fileName);
    %......................................................................
    %Rutina ode15s
    fprintf('\n')
    fprintf(sprintf('Estadísticas, ode15s, Gierer-Meinhardt, tiempo= %d, Rectangulo\n',tf));
    fprintf('------------\n');
    [t,w] = ode15s(@rhs_turing_Gierer_Meinhardt,[t0,tf],w0,options);
    usol = reshape(w(length(t),1:(N+2)^2),N+2,N+2);
    vsol = reshape(w(length(t),(N+2)^2+1:2*(N+2)^2),N+2,N+2);
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %Rectangulo
    %Componente U
    figure;
    contourf(x,y,usol');
    colorbar;    
    daspect([1 3 1]);
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$u(x,y)$','Interpreter','latex');
    title(['Solución numérica u, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
    set(gca,'CLim',[-0.5 3.5]);
    axis([ax bx ay by -0.5 3.5]);
    fileName=sprintf('./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Rectangulo/Rutina-ode15s/Componente-U/Tiempo-%d-turing.png',tf);
    saveas(gcf,fileName);
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %Componente V      
    figure;
    contourf(x,y,vsol');
    colorbar;    
    daspect([1 3 1]);
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$v(x,y)$','Interpreter','latex');
    title(['Solución numérica v, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
    set(gca,'CLim',[-0.5 3.5]);
    axis([ax bx ay by -0.5 3.5]);
    fileName=sprintf('./Evaluaciones-4.2-C1-C2/Gierer-Meinhardt/Rectangulo/Rutina-ode15s/Componente-V/Tiempo-%d-turing.png',tf);
    saveas(gcf,fileName);
    %......................................................................
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %----------------------------------------------------------------------
    %Schnakenberg
    u_star = (a+b);
    v_star = b/(u_star^2);
    u0 = u_star + 10^-3 * exp(-100*((X-(bx-ax)/3).^2+(Y-(by-ay)/3).^2));
    v0 = v_star * ones(N+2);    
    u0 = reshape(u0,(N+2)*(N+2),1);
    v0 = reshape(v0,(N+2)*(N+2),1);
    w0 = [u0;v0];
    %......................................................................
    %Rutina ode45
%     fprintf('\n')
%     fprintf(sprintf('Estadísticas, ode45, Schnakenberg, tiempo= %d, Rectangulo\n',tf));
%     fprintf('------------\n');
%     [t,w] = ode45(@rhs_turing_Schnakenberg,[t0,tf],w0,options);
%     usol = reshape(w(length(t),1:(N+2)^2),N+2,N+2);
%     vsol = reshape(w(length(t),(N+2)^2+1:2*(N+2)^2),N+2,N+2);
%     %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%     %Rectangulo
%     %Componente U
%     figure;
%     contourf(x,y,usol');
%     colorbar;    
%     daspect([1 3 1]);
%     xlabel('$x$','Interpreter','latex');
%     ylabel('$y$','Interpreter','latex');
%     zlabel('$u(x,y)$','Interpreter','latex');
%     title(['Solución numérica u, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
%     set(gca,'CLim',[0 3.5]);
%     axis([ax bx ay by 0 3.5]);
%     fileName=sprintf('./Evaluaciones-4.2-C1-C2/Schnakenberg/Rectangulo/Rutina-ode45/Componente-U/Tiempo-%d-turing.png',tf);
%     saveas(gcf,fileName);
%     %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
%     %Componente V      
%     figure;
%     contourf(x,y,vsol');
%     colorbar;    
%     daspect([1 3 1]);
%     xlabel('$x$','Interpreter','latex');
%     ylabel('$y$','Interpreter','latex');
%     zlabel('$v(x,y)$','Interpreter','latex');
%     title(['Solución numérica v, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
%     set(gca,'CLim',[0 3.5]);
%     axis([ax bx ay by 0 3.5]);
%     fileName=sprintf('./Evaluaciones-4.2-C1-C2/Schnakenberg/Rectangulo/Rutina-ode45/Componente-V/Tiempo-%d-turing.png',tf);
%     saveas(gcf,fileName);
    %......................................................................
    %......................................................................
    %Rutina ode15s
    fprintf('\n')
    fprintf(sprintf('Estadísticas, ode15s, Schnakenberg, tiempo= %d, Rectangulo\n',tf));
    fprintf('------------\n');
    [t,w] = ode15s(@rhs_turing_Schnakenberg,[t0,tf],w0,options);
    usol = reshape(w(length(t),1:(N+2)^2),N+2,N+2);
    vsol = reshape(w(length(t),(N+2)^2+1:2*(N+2)^2),N+2,N+2);
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %Rectangulo
    %Componente U
    figure;
    contourf(x,y,usol');
    colorbar;    
    daspect([1 3 1]);
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$u(x,y)$','Interpreter','latex');
    title(['Solución numérica u, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
    set(gca,'CLim',[0 3]);
    axis([ax bx ay by 0 3]);
    fileName=sprintf('./Evaluaciones-4.2-C1-C2/Schnakenberg/Rectangulo/Rutina-ode15s/Componente-U/Tiempo-%d-turing.png',tf);
    saveas(gcf,fileName);
    %''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''''
    %Componente V      
    figure;
    contourf(x,y,vsol');
    colorbar;    
    daspect([1 3 1]);
    xlabel('$x$','Interpreter','latex');
    ylabel('$y$','Interpreter','latex');
    zlabel('$v(x,y)$','Interpreter','latex');
    title(['Solución numérica v, N = ',num2str(N),', t = ',num2str(t(length(t)),'%6.4f')]);
    set(gca,'CLim',[0 3]);
    axis([ax bx ay by 0 3]);
    fileName=sprintf('./Evaluaciones-4.2-C1-C2/Schnakenberg/Rectangulo/Rutina-ode15s/Componente-V/Tiempo-%d-turing.png',tf);
    saveas(gcf,fileName);
    %......................................................................
end
diary off;
close all;