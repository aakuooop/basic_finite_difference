clc
clear all
close all

lx = 4.00;                                  % length of propagation
dx = 0.1;                                   % grid size
dt = [0.1 0.05] ;                           % step size
m = lx/dx;                                  % No.of intervals

U1 = 1;
U2 = 0;
% Grid Generation

for i = 1:m+1
    x(i) = (i-1)*dx;
end



for k = dt
    c = k/dx;
    % initial condition
    
    for i = 1:m+1
        if x(i)>=0 && x(i)<=2
            u(i) = 1;
        else
            u(i) = 0;
        end
    end
    % boundary condition
    u(1) = U1;
    u(m+1) = U2;
    u_predict(1) = U1;
    u_predict(m+1) = U2;
  
    count = 1;
    t = 0;
    while count>0
        t = t+k
        count = 0;
 % Descretized equation according to MacCormack method  
        % predictor 
       for i = 2:m
           u_predict(i) =  u(i)-0.5*c*(u(i+1)^2-u(i)^2);
       end
       
       % corrector
        unew(1) = U1;
        for i = 2:m
            unew(i)=0.5*(u(i)+u_predict(i)-0.5*c*(u_predict(i)^2-u_predict(i-1)^2));
        end
        unew(m+1) = U2;
        u = unew;
        count = count+1;
        
        if t >= 1.8
            break
        end
        
        
    end
    if k == dt(1)
        plot(x,u,'sr-')
    else
        plot(x,u,'*-b')
    end
    hold on
    
end
legend('c=1','c=0.5','Location','West')
xlabel('x')
ylabel('u')

