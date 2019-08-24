% CODE TO EXECUTE THE DYNAMIC SIMULATION FOR A SINGLE INITIAL CONFIGURATION

clc
close all

global K0
% K0 is a switch: 1- all forces at play, 0- only depletion and excluded
% volume interactions are ON.

% Look at the "initial_conditions_Shenetal.m" file to know possible values
% that can be used for variable- 'intl'
intl=6.4;

imj=1;
if isnumeric(intl)==1
    n=floor(intl);
else
    n=10;
end

% Time of simulation
t=300;

% Initial conditions for the drops: with perturbation
initial_conditions_Shenetal
x0=x0+0.2*rand(n,1);
y0=y0+0.2*rand(n,1);

%%
qq=1;
close all

for K0=1
    
    % Code That Solves For The Position and Velocity
    Other_model
    
    ddrop=zeros(t,1);
    for it=1:t
        cent=[mean(x_drop(:,it)) mean(y_drop(:,it))];
        ddrop(it)=sum(sqrt((x_drop(:,it)-cent(1)).^2+(y_drop(:,it)-cent(2)).^2));
    end
    

    figure(1)
    plot(1:t,ddrop)
    hold all
    fig=gcf;
    fig.Units='normalized';
    fig.OuterPosition=[0 0.5 0.3 0.3];
    
    [xc,yc]=cylinder(1,100);

    snapt=ceil(linspace(1,t,5));
    
    if K0>0
        cc=[0.85 0.33 0.1];
    else
        cc=[0 0.5 0.5];
    end
    
    figure(2)
    for q=snapt
        subplot(2,length(snapt),qq)
        qq=qq+1;
        for k=1:n
            fill(xc(1,:)+x_drop(k,q),yc(1,:)+y_drop(k,q),cc)
            hold on
            axis equal
        end
        title(num2str(q))
        axis([-2+min(x_drop(:,q)) 2+max(x_drop(:,q)) -2+min(y_drop(:,q)) 2+max(y_drop(:,q))])
    end
    fig=gcf;
    fig.Units='normalized';
    fig.OuterPosition=[0 0 1 0.4];
end