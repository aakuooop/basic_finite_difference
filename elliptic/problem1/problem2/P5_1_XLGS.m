tic
clc
clear all
close all

% Length of the rectangular plate
lx = 5;

% Width of the rectangular plate
ly = 5;

% Grid size in horizontal direction
dx = 0.2;

% Grid size in vertical direction
dy = 0.2;

% Number of elements along the length of the plate
m = lx/dx;

% Number of elements along the width of the plate
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

% Line Gauss-Seidel Iteration Method given by (5.16)

for i = 1:m-1
    a(i) = 1;
    b(i) = -2*(1+beta*beta);
    c(i) = 1;
end

% Create triadiagonal matrix 
P = diag(a(1:m-2),-1)+diag(b(1:m-1))+diag(c(1:m-2),1);

errormax = 0.01;

count = (m-1)*(n-1);
iteration = 0;
while count>0
    count = 0;
    prev = psi;
    
    for j = 2:n
        p = 0;
        for i = 2:m
            p = p+1;
       % RHS of equation (5.16)
            if (i==2)
                Q(p) = -beta^2*(psi(i,j+1)+psi(i,j-1))-psi(1,j);
            elseif (i==m)
                Q(p) = -beta^2*(psi(i,j+1)+psi(i,j-1))-psi(m+1,j);
            else
                 Q(p) = -beta^2*(psi(i,j+1)+psi(i,j-1));
            end
         
        end
        R = P\Q';
        p = 0;
        for i = 2:m
            p = p+1;
        psi(i,j) = R(p);
        end
        
    end 
   
    error = 0;
    for j = 2:n
        for i = 2:m
            error = error+abs(psi(i,j)-prev(i,j));
        end
    end
    count = count+1;
    iteration = iteration+1, error
    if (error<errormax)
        break;
    end
end
fid = fopen( 'P5_1_XLGS.txt', 'wt' );
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
    
     figure(2)
    contourf(x,y,psi,25)
    






