clc
clear all
close all

% Given input data


a = 300;                                          % speed of sound at the altitude of 10 km
dx =5.0;                                          % Grid size
dt = [0.0167 0.016666 0.015 0.0075];              % step size
lx = 300;                                         % length of one-dimensional tube
m = lx/dx;                                        % No. of intervals
nstep = 32;                                       % No. of time steps
U1 = 0;
U2 = 0;
% Grid generation

for i = 1:m+1
    x(i) = (i-1)*dx;
end


% initial time

for k = dt
    
    c = a*k/dx;                                   % Courant number 
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
    u(1) = U1;                                         % velocity when x = 0 at every timestep
    u(m+1) = U2;                                       % velocity when x = lx at every time step

    t = 0;
    count = 1;
    while count>0
        count = 0;
        t = t+k
        unew(1)=U1;
        for i = 2:m
      % first upwind differencing technique according to equation (6.16)
            unew(i) = u(i)-c*(u(i)-u(i-1)); 

        end
        unew(m+1)=U2;
        u = unew;
        count = count+1;
        if t>=0.45
            break
        end
    end
   
   if k == dt(1)
       plot(x,u,'*-r')
   elseif k == dt(2)
       plot(x,u,'^-r')
   elseif k == dt(3)
       plot(x,u,'o-r')
   else
       plot(x,u,'d-r')
   end
   hold on
end
legend('Exact','c=0.9996','c=0.90','c=0.45','Location','NorthWest')
xlabel('x')
ylabel('u')

 