tic
clc
clear all
close all

% Length of the rectangular plate
lx = 1;

% Width of the rectangular plate
ly = 2;

% Grid size in horizontal direction
dx = 0.05;

% Grid size in vertical direction
dy = 0.05;

% Number of elements along the length of the plate
m = lx/dx;

% Number of elements along the width of the plate
n = ly/dy;

beta = dx/dy;


% Temperature on the bottom end of the plate (boundary condition)
T1 = 100;

for j = 1:n+1
    for i = 1:m+1
        x(i,j) = (i-1)*dx;
        y(i,j) = (j-1)*dy;
    end
end


% initial condition &&  boundary condition

T = zeros(m+1,n+1);

T(1:m+1,1) = T1;



for i = 1:m-1
    a1(i) = 1;
    b1(i) = -2*(1+beta*beta);
    c1(i) = 1;
end

for i = 1:n-1
    a2(i) = beta^2;
    b2(i) = -2*(1+beta*beta);
    c2(i) = beta^2;
end

% Discretized equation (5.22) and (5.33)
% Form tridiagonal matrix

P1 = diag(a1(1:m-2),-1)+diag(b1(1:m-1))+diag(c1(1:m-2),1);
P2 = diag(a2(1:n-2),-1)+diag(b2(1:n-1))+diag(c2(1:n-2),1);

errormax = 0.01;

count = (m-1)*(n-1);
iteration = 0;
while count>0
    count = 0;
    prev = T;
    
    for j = 2:n
        p = 0;
        for i = 2:m
            p = p+1;
            % RHS of equation (5.22)
            if (i==2)
                Q1(p) = -beta^2*(T(i,j+1)+T(i,j-1))-T(1,j);
            elseif (i==m)
                Q1(p) = -beta^2*(T(i,j+1)+T(i,j-1))-T(m+1,j);
            else
                Q1(p) = -beta^2*(T(i,j+1)+T(i,j-1));
            end
            
        end
        R1 = P1\Q1';
        a = 0;
        for i = 2:m
            a = a+1;
            T(i,j) = R1(a);
        end
        
    end
    
    for i = 2:m
        p = 0;
        for j = 2:n
            p = p+1;
            % RHS of equation (5.23)
            if (j==2)
                Q2(p) = -(T(i-1,j)+T(i+1,j))-T(i,1);
            elseif (j==n)
                Q2(p) = -(T(i-1,j)+T(i+1,j))-T(i,n+1);
            else
                Q2(p) = -(T(i-1,j)+T(i+1,j));
            end
            
        end
        R2 = P2\Q2';
        a = 0;
        for j = 2:n
            a = a+1;
            T(i,j) = R2(a);
        end
        
    end
    
    error = 0;
    for j = 2:n
        for i = 2:m
            error = error+abs(T(i,j)-prev(i,j));
        end
    end
    count = count+1;
    iteration = iteration+1, error
    if (error<errormax)
        break;
    end
end

%contourf(x,y,T)
%Plotting
surf(x,y,T);view(2)
shading interp
grid on
colorbar
title(' Steady state temperature distribution on a two-dimensional rectangular plate')
xlabel('Length of the plate')
ylabel('Width of the plate')
% disp('        Y       X = 0.0   X = 0.1   X = 0.2   X = 0.3   X = 0.4   X = 0.5   X = 0.6   X = 0.7   X = 0.8   X = 0.9   X = 1.0 ');

% Opening new file and writing data on it
fid = fopen( 'ADI_temperature_distribution.txt', 'wt' );
fprintf(fid,'        Y       X = 0.0   X = 0.1   X = 0.2   X = 0.3   X = 0.4   X = 0.5   X = 0.6   X = 0.7   X = 0.8   X = 0.9   X = 1.0 \n\n');
for j = n+1:-1:1
    A = [y(1,j);T(1,j);T(3,j);T(5,j);T(7,j);T(9,j);T(11,j);T(13,j);T(15,j);T(17,j);T(19,j);T(21,j)];
    fprintf(fid,' %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n',A);
end

toc






