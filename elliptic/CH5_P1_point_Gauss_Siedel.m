
clc
clear all
close all

lx = 1;
ly = 2;
dx = 0.05;
dy = 0.05;
m = lx/dx;
n = ly/dy;
T1 = 100;
beta = dx/dy;
niter = 800;        % No. of iterations 

% Grid 
for i = 1:m+1
    for j = 1:n+1
        x(i,j) = (i-1)*dx;
        y(i,j) = (j-1)*dy;
    end
end


% initial condition &&  boundary condition


T = zeros(m+1,n+1);           

T(1:m+1,1) = T1;      % Boundary condition

errormax = 0.01;

% Discretized equation (5.15)

for k = 1:niter
     prev = T;
    for i = 2:m
        for j = 2:n
            T(i,j) = 1/(2*(1+beta*beta))*(T(i+1,j)+T(i-1,j)+beta*beta*(T(i,j+1)+T(i,j-1)));
        end
    end 
    error = 0 ;
    
    % Convergence Criteria given by equation (5.27)
    for i = 2:m
        for j = 2:n
            error = error + abs(T(i,j)-prev(i,j)); 
        end
    end
    k , error 
    if (error<errormax)
        break;
    end
end

% Temperature contour 
% contour(x,y,T,20)
surf(x,y,T);view(2)
shading interp
grid on
colorbar
xlabel('X')
ylabel('Y')

% Opening file and writing data on it 
fid = fopen( 'temperature_distribution.txt', 'wt' );
fprintf(fid,'        Y       X = 0.0   X = 0.1   X = 0.2   X = 0.3   X = 0.4   X = 0.5   X = 0.6   X = 0.7   X = 0.8   X = 0.9   X = 1.0 \n\n');
 for j = n+1:-1:1
      A = [y(1,j);T(1,j);T(3,j);T(5,j);T(7,j);T(9,j);T(11,j);T(13,j);T(15,j);T(17,j);T(19,j);T(21,j)];
      fprintf(fid,' %10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f%10.2f\n',A);
 end 








