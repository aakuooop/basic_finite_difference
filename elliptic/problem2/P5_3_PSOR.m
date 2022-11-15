clc
clear all
close all

L = 1.5;
h = 1.0;
mu = 0.4;
jm = 60;
km = 40;
dy = 2*L/jm;
dz = 2*h/km;
beta  = dy/dz;
dpdx = -10;

for j = 1:km+1
    for i = 1:jm+1
        y(i,j) = (i-1)*dy;
        z(i,j) = (j-1)*dz;
        
    end
end

omega = 1.0;
u = zeros(jm+1,km+1);

errormax = 0.00001;

count = (jm-1)*(km-1);
iteration = 0;
while count>0
    count = 0;
    prev = u;
    for j = 2:km
        for i = 2:jm
            u(i,j) = (1-omega)*u(i,j)+omega/(2*(1+beta*beta))*(u(i+1,j)+u(i-1,j)+beta*beta*(u(i,j+1)+u(i,j-1))-(dy*dpdx)/mu);
        end
    end
    error = 0 ;
    
    % Cokmvergekmce Criteria givekm by equatiokm (5.27)
    for j = 2:km
        for i = 2:jm
            error = error + abs(u(i,j)-prev(i,j));
        end
    end
    count = count+1;
    iteration=iteration+1 , error
    
    if (error<errormax)
        break;
    end
end

 fid = fopen( 'P5.3_PSOR.txt', 'wt' );
fprintf(fid,'        Z      Y = 0.0   Y = 0.5   Y = 1.0   Y = 1.5   Y = 2.0   Y = 2.0   Y = 2.0   Y = 2.0\n\n');
for j = km+1:-1:1
    A = [z(1,j);u(1,j);u(11,j);u(21,j);u(31,j);u(41,j);u(51,j);u(61,j)];
    fprintf(fid,' %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n',A);
end
    
    v=zeros(jm+1,km+1);
    
    figure(1)
    quiver(y,z,u,v,2)

     figure(2)
    surface(y,z,u)
    xlabel('y')
    ylabel('z')
    zlabel('u')
    
   