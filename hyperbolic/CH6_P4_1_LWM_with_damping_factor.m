clc
clear all
close all

% Given input data


a = 300;                        % speed of sound at the altitude of 10 km
dx =5.0;                        % Grid size
dt = 0.015;                     % step size
c = a*dt/dx;                    % Courant number
lx = 300;                       % length of one-dimensional tube
m = lx/dx;                      % No. of intervals
U1 = 0;
U2 = 0;
% Grid generation

for i = 1:m+1
    x(i) = (i-1)*dx;
end

nepsil = 0;

for epsilon = 0.03:0.03:0.12   % various damping factor 
    
    % initial condition
    
    for i = 1:m+1
        
        if x(i)>=0 && x(i)<=50
            u(i) = 0;
        elseif x(i)>=50 && x(i)<=110
            u(i) = 100*sin(pi*((x(i)-50)/60));
        else
            u(i) =  0;
        end
        
    end
    
    % boundary condition
    
    
    u(1) = U1;                % velocity when x = 0 at every timestep
    u(m+1) = U2;              % velocity when x = lx at every timestep
    
    t = 0;                                             % initial time
    nstep = 0;
    
    count = 1;
    nepsil = nepsil+1;
    subplot(2,2,nepsil)
    
    plot(x,u,'r*-')
    hold on
    while count>0
        count = 0;
        t = t+dt;
        unew(1) = U1;
        
      % Discretized equation according to the Lax-Wendroff method with damping factor   
        for i = 2:m
            
            unew(i) = u(i)-c*0.5*(u(i+1)-u(i-1))+0.5*c^2*(u(i+1)-2*u(i)+u(i-1))+epsilon*(u(i+1)-2*u(i)+u(i-1));
            
        end
        unew(m+1) = U2;
        u = unew;
        count = count+1;
        nstep = nstep+1;
        
        
        
        if nstep==12
            plot(x,u,'r+-')
        elseif nstep==24
            plot(x,u,'rd-')
        elseif nstep==36
            plot(x,u,'rs-')
        end
        hold on
        
        if nstep==36
            break
        end
        
        
    end
    legend('T=0.0','T=0.18','T=0.36','T=0.54','Location','East')
    axis ([0 350 -20 120])
    txt = ['\epsilon = ',num2str(epsilon)];
    text(175, -37, txt,'FontSize',10)
    
end




