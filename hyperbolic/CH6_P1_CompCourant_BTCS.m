clc
clear all
close all

% Given input data


a = 300;                                          % speed of sound at the altitude of 10 km
dx =5.0;                                          % Grid size
dt = [0.01666666667 0.015 0.03];                  % step size
lx = 300;                                         % length of one-dimensional tube
m = lx/dx;                                        % No. of intervals
U1 = 0;
U2 = 0;
% Grid generation

for i = 1:m+1
    x(i) = (i-1)*dx;
end

% Fully implicit Euler's BTCS method given by equation (6.6)

for k = dt
    k
    c = a*k/dx
    
    for i = 1:m-1
        d1(i)=0.5*c;
        d2(i)=-1;
        d3(i)=-0.5*c;
    end
    
 % Formulate tridiagonal matrix 
    P = diag(d1(1:m-2),-1)+diag(d2(1:m-1))+diag(d3(1:m-2),1);
    
    for i = 1:m+1
        
        if x(i)>=0 && x(i)<=50
            u(i) = 0;
        elseif x(i)>=50 && x(i)<=110
            u(i) = 100*sin(pi*((x(i)-50)/60));
        else
            u(i) =  0;
        end
        
    end
    
    u(1) = U1;                                         % velocity when x = 0 at all timesteps
    u(m+1) = U2;                                       % velocity when x = lx at all time steps
    
    t = 0;
    count = 1;
    while count>0
        count = 0;
        t = t+k;
        r = 0;
        for i = 2:m
            r = r+1;
            if i ==2
                Q(r) = -u(i)-0.5*c*u(i-1);
            elseif i==m
                Q(r) = -u(i)+0.5*c*u(i+1);
            else
                Q(r) = -u(i);
            end
        end
     % solving matrix
        R = P\Q';
        r = 0;
        for i = 2:m
            r = r+1;
            u(i)= R(r);
        end
        
        count = count+1;
        if t>=0.45
            break
        end
    end
    t
    if k == dt(1)
        plot(x,u,'*-r')
    elseif k == dt(2)
        plot(x,u,'^-r')
    elseif k == dt(3)
        plot(x,u,'o-r')
    end
     
    hold on
end
legend('Exact','c=0.9','c=1.8','Location','NorthWest')
xlabel('x')
ylabel('u')


















