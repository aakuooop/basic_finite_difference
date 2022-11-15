
clc
clear all
close all

lx = 5;
ly = 5;
dx = 0.2;
dy = 0.2;
m = lx/dx;
n = ly/dy;
beta = dx/dy;


for j = 1:n+1
    for i = 1:m+1
        x(i,j) = (i-1)*dx;
        y(i,j) = (j-1)*dy;
    end
end

psi = zeros(m+1,n+1);

for j = 1:n+1
    for i = 1:m+1
        if x(i,j)<=1 && y(i,j)==0
            psi(i,j) = 0.0;
        elseif x(i,j)>=1.2 && y(i,j)==0
            psi(i,j) = 100.0;
        elseif x(i,j)==0 && y(i,j)<=ly
            psi(i,j) = 0.0;
        elseif x(i,j)==lx && y(i,j)<=3
            psi(i,j) = 100.0;
        elseif x(i,j)==lx && y(i,j)>=3.2
            psi(i,j) = 0.0;
        elseif x(i,j)>=0 && y(i,j)==ly
            psi(i,j) = 0.0;
        else
            psi(i,j) = 0.0;
        end
    end
end

omega =  1.785;

iteration = 0;
count = (m-1)*(n-1);
errormax = 0.01;
while count>0
    count = 0;
    prev = psi;
  
    % Discritised equation (5.18)
    for i = 2:m
        for j = 2:n
            psi(i,j) = (1-omega)*psi(i,j)+omega/(2*(1+beta*beta))*(psi(i+1,j)+psi(i-1,j)+beta*beta*(psi(i,j+1)+psi(i,j-1)));
        end
    end 
   
   error = 0;
   
   % Convergence criterion 
    for i = 2:m
        for j = 2:n
            error = error+abs(psi(i,j)-prev(i,j));
        end
    end
    count = count+1;
    iteration = iteration+1 , error
    if (error<errormax)
        break;
    end
end

fid = fopen( 'P5_1_PSOR.txt', 'wt' );
fprintf(fid,'        Y       X = 0.0   X = 1.0   X = 2.0   X = 3.0   X = 4.0   X = 5.0 \n\n');
for j = n+1:-1:1
    A = [y(1,j);psi(1,j);psi(6,j);psi(11,j);psi(16,j);psi(21,j);psi(26,j)];
    fprintf(fid,' %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n',A);
end

u = zeros(n+1,m+1);
v = zeros(n+1,m+1);

for j = 2:n
    for i = 2:m
           
           u(i,j) = -(psi(i,j+1)-psi(i,j-1))/(2*dy);
            v(i,j) = (psi(i+1,j)-psi(i-1,j))/(2*dx);
      
    end
end

    figure(1)
     quiver(x,y,u,v,2)
     hold on 
    contour(x,y,psi,25)







