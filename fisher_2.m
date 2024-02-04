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
% rmdir("./Evaluaciones");
mkdir("Evaluaciones");
%__________________________________________________________________________
%Valores para distintas evaluaciones
tiempos=[0.5 1 2 10];
%__________________________________________________________________________
%DATOS DEL PROBLEMA

global h N
a = -2 ;
b = 2 ;

N = 47;             % Número de nodos interiores en las direcciones x e y
h = (b-a)/(N+1) ;             % Tamaño de malla
 
x = linspace(a,b,N+2);
y = linspace(a,b,N+2);

t0 = 0 ;
tf = 10 ;
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
%FORs de tipos
for tipo=1:5 %Tipo de difusión
    %__________________________________________________________________________
    
    %--------------------------------------------------------------------------------------------
%     % Representación gráfica del tamaño de paso en tiempo: 
%     figure(1)
%     semilogy(t1(2:end),diff(t1),'o','Markersize',5);
%     xlabel('$t$','Interpreter','latex')
%     ylabel('$\Delta t$','Interpreter','latex')
%     title('Tamaño de paso en tiempo')
%     fileName=sprintf('./Evaluaciones/Rutina-ode45-tipo-difusion-%d-fisher.png',tipo);
%     saveas(gcf,fileName)
%     %--------------------------------------------------------------------------------------------
    %For tiempos
    for temp = 1:4
        i=tiempos(temp);

%         %Rutina ode45
%         fprintf('\n')
%         fprintf('Estadísticas\n');
%         fprintf('------------\n');
%         [t1,u1] = ode45(@(t,u)rhs_fisher2(t,tipo,u),[t0,i],u0,options); %ode45
% 
%         
%         fileName=sprintf('./Evaluaciones/Rutina-ode45-tiempo-%d-tipo-difusion-%d-fisher.png',temp,tipo);
%         sol = reshape(u1(length(t1),:),N+2,N+2);
%         figure(1)
%         surf(x,y,sol');
%         %contourf(x,y,sol');
%         xlabel('$x$','Interpreter','latex');
%         ylabel('$y$','Interpreter','latex');
%         zlabel('$u(x,y)$','Interpreter','latex');
%         title(['Solución numérica, N = ',num2str(N),', t = ',num2str(i,'%6.4f')]);
%         colorbar;
%         colormap('jet');
%         %set(gca,'CLim',[0 1.3]);
%         axis([a b a b 0 1.3]);
%         saveas(gcf,fileName)

        %__________________________________________________________________________
        %Rutina ode15sfprintf('\n')
        fprintf('Estadísticas\n');
        fprintf('------------\n');
        [t2,u2] = ode15s(@(t,u)rhs_fisher2(t,tipo,u),[t0,i],u0,options); %ode15s
        %--------------------------------------------------------------------------------------------

        fileName=sprintf('./Evaluaciones/Rutina-ode15s-tiempo-%d-tipo-difusion-%d-fisher.png',temp,tipo);
        sol = reshape(u2(length(t2),:),N+2,N+2);
        figure(2)
        surf(x,y,sol');
        %contourf(x,y,sol');
        xlabel('$x$','Interpreter','latex');
        ylabel('$y$','Interpreter','latex');
        zlabel('$u(x,y)$','Interpreter','latex');
        title(['Solución numérica, N = ',num2str(N),', t = ',num2str(i,'%6.4f')]);
        colorbar;
        colormap('jet');
        set(gca,'CLim',[0 1.3]);
        axis([a b a b 0 1.3]);
        saveas(gcf,fileName)
    end
    % Representación gráfica del tamaño de paso en tiempo ode45: 
%     figure(3)
%     semilogy(t1(2:end),diff(t1),'o','Markersize',5);
%     xlabel('$t$','Interpreter','latex')
%     ylabel('$\Delta t$','Interpreter','latex')
%     title('Tamaño de paso en tiempo')
%     fileName=sprintf('./Evaluaciones/Rutina-ode45-tipo-difusion-%d-fisher.png',tipo);
%     saveas(gcf,fileName)
    
%     % Representación gráfica del tamaño de paso en tiempo ode15s: 
%     figure(4)
%     semilogy(t2(2:end),diff(t2),'o','Markersize',5);
%     xlabel('$t$','Interpreter','latex')
%     ylabel('$\Delta t$','Interpreter','latex')
%     title('Tamaño de paso en tiempo')
%     fileName=sprintf('./Evaluaciones/Rutina-ode15s-tipo-difusion-%d-fisher.png',tipo);
%     saveas(gcf,fileName)
%     %--------------------------------------------------------------------------------------------    
end