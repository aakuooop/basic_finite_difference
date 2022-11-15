clc
clear all
close all

% Given input data


a = 300;                                          % speed of sound at the altitude of 10 km
dx =5.0;                                          % Grid size
dt = 0.015;                                       % step size
c = a*dt/dx;                                      % Courant number
lx = 300;                                         % length of one-dimensional tube
m = lx/dx;                                        % No. of intervals
U1 = 0;
U2 = 0;
ep1 = 0.026;
ep2 = 0.002;
% Grid generation

for i = 1:m+1
    x(i) = (i-1)*dx;
end


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


u(1) = U1;                                         % velocity when x = 0 at all timesteps
u(m+1) = U2;                                       % velocity when x = lx at all time steps

u_predict(1) = 1;
u_predict(m+1) = 0;



t = 0;                                             % initial time
nstep = 0;
count = 1;
plot(x,u,'r*-')
hold on
while count>0
    count = 0;
    t = t+dt
    u_correct(1) = U1;
    % Second-order accurate Lax-Wendroff method with flux corrected transport
%given by equation (6.38a) and (6.38b)
     % predictor 
       for i = 2:m
           u_predict(i) =  u(i)-0.5*c*(u(i+1)-u(i-1))+(ep1+0.5*c^2)*(u(i+1)-2*u(i)+u(i-1));
       end
    % corrector
        for i = 2:m
           u_correct(i) = u_predict(i)-ep2*(u_predict(i+1)-2*u_predict(i)+u_predict(i-1));
        end
    u_correct(m+1) = U2;
    u = u_correct;
    count = count+1;
    nstep = nstep+1;
    if nstep==12
        plot(x,u,'rd-')
    elseif nstep==24
        plot(x,u,'rs-')
    elseif nstep==36
        plot(x,u,'r-')
    end
    hold on
    
    if nstep==36
        break
    end
    
    
end


legend('T=0.0','T=0.18','T=0.36','T=0.54','Location','East')
axis ([0 350 -20 120])















