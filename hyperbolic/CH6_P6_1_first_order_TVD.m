clc
clear all
close all

lx = 4.00;                                  % length of propagation
dx = 0.1;                                   % grid size
dt = 0.1;                                   % step size
m = lx/dx;                                  % No.of intervals
nstep = 18;                                 % Number of time step
c = dt/dx
U1 = 1;
U2 = 0;
% Grid Generation

for i = 1:m+1
    x(i) = (i-1)*dx;
end


% initial condition

for i = 1:m+1
    if x(i)>=0 && x(i)<=2
        u(1,i) = 1;
    else
        u(1,i) = 0;
    end
end


% boundary condition
for k = 1:nstep+1
    u(k,1) = U1;
    u(k,m+1) = U2;
end

u_predict(1) = U1;
u_predict(m+1) = U2;


t = 0;


for k = 1:nstep
    t = t+dt
   
    for i = 1:m+1
        E(i) = 0.5*u(k,i)^2;
    end
 % Descretized equation according to first order TVD scheme given by
 % equation (6.56a)and (6.56b) 
    
    % predictor
    for i = 2:m
        u_predict(i) =  u(k,i)-0.5*c*(E(i+1)-E(i-1));  %(6.56a)
    end
    % corrector
    
    
    for i = 2:m
        delu_for = u(k,i+1)-u(k,i);
        delu_bac = u(k,i)-u(k,i-1);
        if delu_for == 0
            alpha_for = u(k,i);
        else 
            alpha_for = (E(i+1)-E(i))/(u(k,i+1)-u(k,i));
        end
        
        if delu_bac == 0
            alpha_bac = u(k,i);
        else
            alpha_bac = (E(i)-E(i-1))/(u(k,i)-u(k,i-1));
        end
        phi_for = abs(alpha_for)*delu_for;
        phi_bac = abs(alpha_bac)*delu_bac;
        
        u(k+1,i) = u_predict(i)+0.5*c*(phi_for-phi_bac); %(6.56b)
    end
end

fid = fopen('CH6_P6_1_First order TVD_solution.txt','wt');
fprintf(fid,'       x      t=0.0     t=0.3     t=0.6     t=0.9     t=1.2     t=1.5     t=1.8\n\n');
for i = 1:m+1
    A = [x(i);u(1,i);u(4,i);u(7,i);u(10,i);u(13,i);u(16,i);u(19,i)];
    fprintf(fid,'%10.2f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f%10.5f\n',A);
end

plot(x,u(1,:),'r.-')
hold on
plot(x,u(4,:),'r.-')
hold on
plot(x,u(7,:),'r.-')
hold on
plot(x,u(10,:),'r.-')
hold on
plot(x,u(13,:),'r.-')
hold on
plot(x,u(16,:),'r.-')
hold on
plot(x,u(19,:),'r.-')


legend('T=0.0','T=0.3','T=0.6','T=0.9','T=1.2','T=1.5','T=1.8','Location','West')
