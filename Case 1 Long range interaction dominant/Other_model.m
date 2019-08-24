% Code That Solves For The Position and Velocity

global n
parameters_final

% Reduce the time step incase of numerical instability
t_step=0.5;

%initialization
x_drop=zeros(n,t+1);
y_drop=zeros(n,t+1);
 
% intl='rand1';

x_drop(:,1)=x0;
y_drop(:,1)=y0;

%--------------------------------------------------------------------
% Flow at infinity
e_infx=1;
e_infy=0;
e_inf=[e_infx,e_infy];

for i=1:t

    r0 = [x_drop(:,i)' y_drop(:,i)'];
    
    %POSITION UPDATE
    [T, P] = ode45(@position,[(i-1)*t_step i*t_step],r0);
    x_drop(:,i+1)=P(end,1:n);
    y_drop(:,i+1)=P(end,n+1:2*n);
    
end
