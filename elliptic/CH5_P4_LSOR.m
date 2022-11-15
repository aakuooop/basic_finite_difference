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
errormax = 0.01;
fid = fopen('LSOR_relaxation_table.txt','wt');
fprintf(fid,' Relaxation parameter    Number of iterations\n\n');

% omega is a relaxation parameter 
for omega = 1.00:0.01:1.32
    omega
    % initial condition &&  boundary condition
    T = zeros(m+1,n+1);
    T(1:m+1,1) = T1;
    % Form tridiagonal matrix
    for i = 1:m-1
        a(i) = omega;
        b(i) = -2*(1+beta*beta);
        c(i) = omega;
    end
    P = diag(a(1:m-2),-1)+diag(b(1:m-1))+diag(c(1:m-2),1);
    count = (m-1)*(n-1);
    iteration = 0;
    while count>0
        count = 0;
        prev = T;
        
        for j = 2:n
            p = 0;
            for i = 2:m
                p = p+1;
                % RHS of equation (5.21)
                if (i==2)
                    Q(p) = -(1-omega)*2*(1+beta^2)*T(i,j)-omega*beta^2*(T(i,j+1)+T(i,j-1))-T(1,j);
                elseif (i==m)
                    Q(p) = -(1-omega)*2*(1+beta^2)*T(i,j)-omega*beta^2*(T(i,j+1)+T(i,j-1))-T(m+1,j);
                else
                    Q(p) = -(1-omega)*2*(1+beta^2)*T(i,j)-omega*beta^2*(T(i,j+1)+T(i,j-1));
                end
                
            end
            R = P\Q';
            p = 0;
            for i = 2:m
                p = p+1;
                T(i,j) = R(p);
            end
            
        end
        
        error = 0;
        for j = 2:n
            for i = 2:m
                error = error+abs(T(i,j)-prev(i,j));
            end
        end
        count = count+1;
        iteration = iteration+1; error;
        if (error<errormax)
            break;
        end
    end
    iteration
    fprintf(fid,'%10.2f\t\t\t\t\t%10.2f\n',omega,iteration);
    plot(iteration,omega,'*-')
    hold on
    
end
xlabel('Number of iteration')
 ylabel('\omega')

toc






