tic
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


for i = 1:m+1
    for j = 1:n+1
        x(i,j) = (i-1)*dx;
        y(i,j) = (j-1)*dy;
    end
end


% initial condition &&  boundary condition

T = zeros(m+1,n+1);

T(1:m+1,1) = T1;
T_init = T;
errormax = 0.01;
fid = fopen('relaxation_table.txt','wt');
fprintf(fid,' Relaxation parameter    Number of iterations\n\n')

for omega = 0.8:0.01:1.98
    omega
iteration = 0;
count = (m-1)*(n-1);
T = T_init;
error = 1;
while error > errormax
    count = 0;
    prev = T;
  
    % Discritised equation (5.18)
    for i = 2:m
        for j = 2:n
            T(i,j) = (1-omega)*T(i,j)+omega/(2*(1+beta*beta))*(T(i+1,j)+T(i-1,j)+beta*beta*(T(i,j+1)+T(i,j-1)));
        end
    end 
   
   error = 0;
   
   % Convergence criterion 
    for i = 2:m
        for j = 2:n
            error = error+abs(T(i,j)-prev(i,j));
        end
    end
    
    iteration = iteration+1;  error;

end
iteration
if omega>=1.00
fprintf(fid,'%10.2f\t\t\t\t\t%10.2f\n',omega,iteration);
end
plot(iteration,omega,'*')
hold on

end

 xlabel('Number of iteration')
 ylabel('\omega')


toc






