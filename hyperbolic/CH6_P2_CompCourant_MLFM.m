clc
clear all
close all

% Given input data


a = 300;                                          % speed of sound at the altitude of 10 km
dx =5.0;                                          % Grid size
dt = [0.01666 0.015 0.0075];                      % step size
lx = 300;                                         % length of one-dimensional tube
m = lx/dx;                                        % No. of intervals
U1 = 0;
U2 = 0;
maxstep = 65;
% Grid generation

for i = 1:m+1
    x(i) = (i-1)*dx;
end



for k = 1:maxstep+1
    u(k,1) = 0;                                         % velocity when x = 0 at all timesteps
    u(k,m+1) = 0;                                       % velocity when x = lx at all time steps
end

for k = dt
    
    c = a*k/dx;
    % initial condition
    
    for i = 1:m+1
        
       if x(i)>=0 && x(i)<=100
            u(1,i) = 0;
        elseif x(i)>=50 && x(i)<=220
            u(1,i) = 100*sin(pi*((x(i)-100)/120));
        else
            u(1,i) =  0;
        end
        
    end
    
    
    
    t = 0;    % initial time 
    
    for  j= 1:maxstep+1
        
        t = t+k
   % Midpoint leapfrog method given by equation (6.18)  
        for i = 2:m
            if j == 1
                u(j+1,i)=u(j,i)+0.5*c^2*(u(j,i-1)-2*u(j,i)+u(j,i+1));
            else
                u(j+1,i)=2*u(j,i)-u(j-1,i)+c^2*(u(j,i-1)-2*u(j,i)+u(j,i+1));
            end
        end
        j
        if t>=0.28
            break
        end
    end
    
    if k == dt(1)
        plot(x,u(j+1,:),'*-r')
    elseif k == dt(2)
        plot(x,u(j+1,:),'^-r')
    else
        plot(x,u(j+1,:),'o-r')
    end
    hold on
end
legend('c=0.9996','c=0.90','c=0.45','Location','North')
xlabel('x')
ylabel('u')

















