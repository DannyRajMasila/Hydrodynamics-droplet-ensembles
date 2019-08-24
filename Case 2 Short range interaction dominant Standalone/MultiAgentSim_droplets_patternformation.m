% CODE FOR SIMULATING DROPLETS IN CHANNELS BASED ON MODEL FROM 
% DANNY AND RENGASWAMY, microfluid nanofluid, 2014

clear
close all

clc
parameters_pf

R1=1.25e-4; % Half of inlet channel dimension
x0=0; %Initial position of drop1 wrt 'x'
y0=0; %Initial position of drop1 wrt 'y'

%INITIAL SPACING
L0=6;
Q1=400;
Q=40;

time=0.25; % SIMULATION TIME

% GA optimized conditions for results reported in the MS: SELECT THE APPROPRIATE L1
L1=[3.0000; 3.7835; 4.1577; 3.9874; 8.1929]; % bad circle
% L1=[5.2846; 3.8021; 4.0756; 3.1889; 4.7573; 7.4025]; % circle with 1 preamble
% L1=[5.4567; 6.6302; 3.9759; 3.5844; 4.2430; 3.0061]; % triangle with 1 preamble
% L1=[4.5844; 4.6386; 4.0018; 3.000; 4.4322; 4.5643; 6.1098; 8.5273]; % square with 1 preamble

n=length(L1)+1;

% Reduce the time step incase of numerical instability
t_step= 1e-4; % single time step

t=round(time/t_step);
% Need for a smaller time step than Danny's code

vx_guess=0.11;
vy_guess=0.01;
y_ec_u=1*R1;
y_ec_l=-1*R1;

v_max=max(1e-9*(Q1+Q)/((y_ec_u-y_ec_l)*(2*R1))/60);


%--------------------------------------------------------------------

drop_in_test_initial=zeros(t,length(L0));
drop_in_test_final=zeros(t,length(L0));
vxd=zeros(n,t,length(L0));
vyd=zeros(n,t,length(L0));
xd=zeros(n,t+1,length(L0));
yd=zeros(n,t+1,length(L0));
t_sim=zeros(length(L0),1);
t_sim_steps=zeros(length(L0),1);


%CHANNEL GEOMETRY
y_ec_u=1*R1;
y_ec_l=-1*R1;

slope_wall=0;

x_b=linspace(0,W,1002);
y_b_u=(heaviside(x_b+R1)-heaviside(x_b-alp*W)).*((h*alp*W-h*x_b+H*x_b)/(alp*W))+(heaviside(x_b-alp*W)-heaviside(x_b-1.2*W)).*(H+slope_wall*(x_b-alp*W));
% y_b_u=(heaviside(x_b+R1)-heaviside(x_b-alp*W)).*(-alp*W*H*h./((H-h).*x_b-alp*W*H))+(heaviside(x_b-alp*W)-heaviside(x_b-1.2*W)).*(H+slope_wall*(x_b-alp*W));

y_b_l=-y_b_u;
slope_b_u=[ones(floor(length(x_b)/2),1);-ones(ceil(length(x_b)/2),1)];

% ADDING ROUGHNESS TO WALL
% Random perturbations with asymmetric disturbance to break top down
% symmetry
y_b_u=y_b_u+0.001*R1*rand(1,length(y_b_u)); 
y_b_l=y_b_l-0.0001*R1*rand(1,length(y_b_l));

y_b_u(1)=y_ec_u;
y_b_l(1)=y_ec_l;
y_b_u(length(y_b_u))=y_ec_u;
y_b_l(length(y_b_l))=y_ec_l;
slope_b_l=-slope_b_u;

lim=ones(t,length(L0));

v_solution = zeros(n,3);

for p=1:length(L0)
    
    Q2=Q(p);
    flowrate=((Q1+Q2)*1e-6*1e-3)./(60*(2*R1));
    
    %initialization
    vx_drop=zeros(n,t);
    vy_drop=zeros(n,t);
    vd_drop=zeros(n,t);
    x_drop=zeros(n,t+1);
    y_drop=zeros(n,t+1);
    sim_time=zeros(t,1);
    
    for k=2:n
        x_drop(k,1)=x_drop(k-1,1)-L1(k-1)*R;
        y_drop(k,1)=y0;
    end
    
    %INITIAL VELOCITIES-GUESS
    v_guess=[vx_guess*ones(n,1);vy_guess*ones(n,1)];
    
    %VELOCITY IN THE INLET AND OUTLET CHANNELS
    v_inletoutlet=1e-9*(Q1+Q2)/((y_ec_u-y_ec_l)*(2*R1))/60;
    
    % %QUASI STEADY STATE SIMULATION------------------------------------------
    drop_in_test=zeros(t,1);
    % name1=strcat('Bifurcation_L0_',num2str(L0(p)),'.mat');
    rer=zeros(100,t);
    
    tic
    
    for i=1:t
        
        
        counter=0;
        counter1=0;
        
        drops_test_temp=zeros(n,1);
        drops_nontest_temp=zeros(n,1);
        
        for u=1:n
            if x_drop(u,i)>=min(x_b)-3*R && x_drop(u,i)<=max(x_b)+3*R
                counter=counter+1;
                drops_test_temp(counter)=u;
            else
                counter1=counter1+1;
                drops_nontest_temp(counter1)=u;
            end
        end
        drops_test=drops_test_temp(1:counter);
        drops_nontest=drops_nontest_temp(1:counter1);
        
        drop_in_test(i)=counter;
        
        drop_in_test_initial(i,p)=drops_test(1);
        drop_in_test_final(i,p)=drops_test(counter);
        
        %VELOCITIES FOR NON TEST SECTION DROPS
        vx_drop(drops_nontest,i)=v_inletoutlet;
        
        start = drops_test(1);
        % only for the test section drops
        if counter>0
            %RE-DOING THE GUESS VALUES FOR THE TEST SECTION DROPS ONLY
            guess=ones(2*counter,1);
            %        for ee=1:counter
            guess(1:counter)=v_guess(drops_test);
            guess(counter+1:2*counter)=v_guess(n+drops_test);
            %        end
            
            % IDTM_non_linear
            % drop_test- contains the identity of drops in the test section
            % counter- number of drops in the test section of the channel
            
            %         parameters
            
            no=zeros(counter,1);
            ino=zeros(10,counter,3);
            inb=zeros(counter,2);
            
            v_sol_big = zeros(5*counter,5*counter);
            g_sol_big = zeros(5*counter,2*counter);
            c_sol_big = zeros(5*counter,1);
            
            y_f_d_big = zeros(5*counter,1);
            d_sol_big = zeros(2*counter,5*counter);
            
            for jc=1:counter
                % identifying drops above and below every drop in the system
                j=drops_test(jc);
                m=x_drop(j,i);
                for o=1:counter
                    oj=drops_test(o);
                    y1=y_drop(oj,i)+sqrt(R^2-(m-x_drop(oj,i))^2);
                    y2=y_drop(oj,i)-sqrt(R^2-(m-x_drop(oj,i))^2);
                    if isreal(y1)==1
                        no(jc)=no(jc)+1;
                        ino(no(jc),jc,1)=oj;
                        ino(no(jc),jc,2)=y1;
                        ino(no(jc),jc,3)=y2;
                    end
                end
                
                %                 ino(no(jc),jc,1)
                
                % identifying boundaries above and below every drop in the system
                
                if m>=min(x_b) && m<=max(x_b)
                    inb(jc,1)=interp1(x_b,y_b_u,m);
                    inb(jc,2)=interp1(x_b,y_b_l,m);
                else
                    inb(jc,1)=y_ec_u;
                    inb(jc,2)=y_ec_l;
                end
            end
            
            rel_err=1;
            err_rate=1;
            ii=1;
            
            while (rel_err>0.5*1e-4) && (abs(err_rate)>0.2*1e-4)
                
                
                
                vx_drop(drops_test,i)=guess(1:counter);
                vy_drop(drops_test,i)=guess(counter+1:2*counter);
                
                vd_drop(:,i)= sqrt(vx_drop(:,i).^2+vy_drop(:,i).^2);
                
                D=zeros(2*counter,2*counter);
                f=zeros(2*counter,1);
                f_d = zeros(2*counter,1);
                
                c_old_1 =0;
                
                for jc=1:counter
                    j=drops_test(jc);
                    m=x_drop(j,i);
                    c=2*no(jc)+2;
                    
                    y_f=zeros(c,2);
                    
                    
                    % No slip conditions
                    y_f(1,1)=inb(jc,1);
                    y_f(1,2)=0;
                    y_f(c,1)=inb(jc,2);
                    y_f(c,2)=0;
                    
                    
                    % Continuity of velocity at every neighbouring drop
                    % (including the drop under consideration)
                    y_f(2:c-1,1)=[ino(1:no(jc),jc,2) ; ino(1:no(jc),jc,3)];
                    y_f(2:c-1,2)=zeros(c-2,1);
                    
                    g_sol_small=zeros(c+1,2*counter);
                    
                    % Linearization
                    g_sol_small(2:c-1,(ino(1:no(jc),jc,1))-start+1)=-beta_value.*[vx_drop(ino(1:no(jc),jc,1),i)/vd_drop(ino(1:no(jc),jc,1),i);...
                        vx_drop(ino(1:no(jc),jc,1),i)/vd_drop(ino(1:no(jc),jc,1),i)];
                    
                    g_sol_small(2:c-1,counter+(ino(1:no(jc),jc,1))-start+1)= -beta_value.*[vy_drop(ino(1:no(jc),jc,1),i)/vd_drop(ino(1:no(jc),jc,1),i);...
                        vy_drop(ino(1:no(jc),jc,1),i)/vd_drop(ino(1:no(jc),jc,1),i)];
                    
                    [a, ix]=sort(y_f(:,1),'descend');
                    temp=zeros(c,1);
                    temp_2 = zeros(c+1,2*counter);
                    
                    % ix stores the order of the eqns. Hence the RHS is
                    % arranged accordingly.
                    for q=1:c
                        temp(q)=y_f(ix(q),2);
                        temp_2(q,:)= g_sol_small(ix(q),:);
                    end
                    y_f(:,1)=a;
                    
                    y_f(:,2)=temp;
                    g_sol_small = temp_2;
                    
                    %                     v4 = ino(1:no(jc),jc,1);
                    
                    % Linearization
                    for v5 = 1:no(jc)
                        g_sol_small(c+1,ino(v5,jc,1)-start+1)= (y_f(2*v5,1)- y_f(2*v5+1,1))*(vx_drop(ino(v5,jc,1),i)/vd_drop(ino(v5,jc,1),i));
                        g_sol_small(c+1,ino(v5,jc,1)+counter-start+1)= (y_f(2*v5,1)- y_f(2*v5+1,1))*(vy_drop(ino(v5,jc,1),i)/vd_drop(ino(v5,jc,1),i));
                    end
                    
                    %Constructing the Flow solution Matrix (v_sol); Size=(c+1)x(c+1)
                    v_sol_small=zeros(c+1,c+1);
                    
                    %'y and 1' terms in c of c+1 rows
                    for g=1:c
                        hh=2;
                        gg=ceil(g/2)*hh;
                        for h=1:2
                            v_sol_small(g,gg)=y_f(g,1)^(2-h);
                            gg=gg+1;
                        end
                    end
                    
                    %First column
                    v_sol_small(1,1)=y_f(1,1)^2;
                    v_sol_small(c,1)=y_f(c,1)^2;
                    v_sol_small(1:c,1)=y_f(:,1).^2;
                    
                    %last row
                    for h=2:2:c
                        v_sol_small(c+1,h)=(y_f(h-1,1)^2-y_f(h,1)^2)/2;
                    end
                    for h=3:2:c+1
                        v_sol_small(c+1,h)=y_f(h-2,1)-y_f(h-1,1);
                    end
                    for h=1:2:c
                        v_sol_small(c+1,1)=v_sol_small(c+1,1)+(y_f(h,1)^3-y_f(h+1,1)^3)/3;
                    end
                    
                    
                    c_sol_small=[y_f(:,2);flowrate];
                    
                    
                    v11= c_old_1 +1;
                    
                    g_sol_big(v11:v11+c,:)= g_sol_small;
                    v_sol_big(v11:v11+c,v11:v11+c)= v_sol_small; % Allocate a big memory for v_sol_big
                    c_sol_big(v11:v11+c,:)=c_sol_small;
                    
                    y_f_d_big(v11:v11+c,1)=[y_f(:,1);0];
                    %----------------------------------------------------------
                    
                    %AVERAGE VELOCITIES
                    p1=find(y_f(:,1)==y_drop(j,i)+R);
                    p2=find(y_f(:,1)==y_drop(j,i)-R);
                    
                    p1=p1(1);
                    p2=p2(1);
                    
                    v1_terms=[1/3*((y_f(p1-1,1)^3)-(y_f(p1,1)^3)), 1/2*((y_f(p1-1,1)^2)-(y_f(p1,1)^2)), 1*(y_f(p1-1,1)-y_f(p1,1))]./(y_f(p1-1,1)-y_f(p1,1));
                    v2_terms=[1/3*((y_f(p2,1)^3)-(y_f(p2+1,1)^3)), 1/2*((y_f(p2,1)^2)-(y_f(p2+1,1)^2)), 1*(y_f(p2,1)-y_f(p2+1,1))]./(y_f(p2,1)-y_f(p2+1,1));
                    
                    % Making of the f matrix [DU=f]
                    f(jc)=0;
                    f(counter+jc)=0;
                    
                    A_term_x = kf/2*(v1_terms(1,1)+v2_terms(1,1));
                    B1_term_x = kf/2*v1_terms(1,2);
                    C1_term_x = kf/2*v1_terms(1,3);
                    B2_term_x = kf/2*v2_terms(1,2);
                    C2_term_x = kf/2*v2_terms(1,3);
                    
                    A_term_y = kf*k1*(v1_terms(1,1)-v2_terms(1,1));
                    B1_term_y = kf*k1*v1_terms(1,2);
                    C1_term_y = kf*k1*v1_terms(1,3);
                    B2_term_y = -kf*k1*v2_terms(1,2);
                    C2_term_y = -kf*k1*v2_terms(1,3);
                    
                    d_sol_big(jc,c_old_1+1)= A_term_x;
                    d_sol_big(jc,c_old_1+p1:c_old_1+p1+1)=[B1_term_x, C1_term_x];
                    d_sol_big(jc,c_old_1+p2+1:c_old_1+p2+2)=[B2_term_x, C2_term_x];
                    
                    d_sol_big(jc+counter,c_old_1+1)= A_term_y;
                    d_sol_big(jc+counter,c_old_1+p1:c_old_1+p1+1)=[B1_term_y, C1_term_y];
                    d_sol_big(jc+counter,c_old_1+p2+1:c_old_1+p2+2)=[B2_term_y, C2_term_y];
                    %                     d_sol_big(jc+counter,5*(jc-1)+1:5*(jc-1)+5)=[A_term_y, B1_term_y, C1_term_y, B2_term_y, C2_term_y];
                    
                    c_old_1 = c_old_1+c+1;
                    
                    % Assembling the D matrix [DU=f]
                    D(jc,jc)=-beta_value*kf;
                    D(counter+jc,counter+jc)=-beta_value*kf;
                    
                    % Effect of drop-drop interactions incorporated in D [DU=f]
                    
                    for jl=1:counter
                        l=drops_test(jl);
                        d=sqrt((x_drop(j,i)-x_drop(l,i))^2+(y_drop(j,i)-y_drop(l,i))^2);
                        
                        v_approach=([vx_drop(j,i) vy_drop(j,i)]*[x_drop(l,i)-x_drop(j,i);  y_drop(l,i)-y_drop(j,i)]...
                            +[vx_drop(l,i) vy_drop(l,i)]*[-x_drop(l,i)+x_drop(j,i) ; -y_drop(l,i)+y_drop(j,i)]);
                        if d==0 || d>=scale_d || v_approach<0   % APPROACH VELOCITY POSITIVE FOR INTERACTION TO OCCUR- NEEDS VERIFICATION
                        else
                            %                     X component equations
                            D(jc,jc)=D(jc,jc)+((x_drop(l,i)-x_drop(j,i))/d)*(kd/(d-2*R))*(((x_drop(j,i)-x_drop(l,i)))/d);
                            D(jc,jl)=D(jc,jl)+((x_drop(j,i)-x_drop(l,i))/d)*(kd/(d-2*R))*(((x_drop(j,i)-x_drop(l,i)))/d);
                            D(jc,counter+jc)=D(jc,counter+jc)+((y_drop(l,i)-y_drop(j,i))/d)*(kd/(d-2*R))*(((x_drop(j,i)-x_drop(l,i)))/d);
                            D(jc,counter+jl)=D(jc,counter+jl)+((y_drop(j,i)-y_drop(l,i))/d)*(kd/(d-2*R))*(((x_drop(j,i)-x_drop(l,i)))/d);
                            %                     Y component equations
                            D(counter+jc,jc)=D(counter+jc,jc)+((x_drop(l,i)-x_drop(j,i))/d)*(kd/(d-2*R))*(((y_drop(j,i)-y_drop(l,i)))/d);
                            D(counter+jc,jl)=D(counter+jc,jl)+((x_drop(j,i)-x_drop(l,i))/d)*(kd/(d-2*R))*(((y_drop(j,i)-y_drop(l,i)))/d);
                            D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+((y_drop(l,i)-y_drop(j,i))/d)*(kd/(d-2*R))*(((y_drop(j,i)-y_drop(l,i)))/d);
                            D(counter+jc,counter+jl)=D(counter+jc,counter+jl)+((y_drop(j,i)-y_drop(l,i))/d)*(kd/(d-2*R))*(((y_drop(j,i)-y_drop(l,i)))/d);
                        end
                    end
                    
                    % Effect of drop-boundary interactions incorporated in D [DU=f]
                    
                    if x_drop(j,i)>=min(x_b) && x_drop(j,i)<=max(x_b)
                        % interaction for the drops in the 2D section
                        
                        [d1, z1]=min((x_b-x_drop(j,i)).^2+(y_b_u-y_drop(j,i)).^2);
                        [d2, z2]=min((x_b-x_drop(j,i)).^2+(y_b_l-y_drop(j,i)).^2);
                        d1=sqrt(d1);
                        d2=sqrt(d2);
                        v_approach_1=[vx_drop(j,i) vy_drop(j,i)]*[(-x_drop(j,i)+x_b(z1))/d1; (-y_drop(j,i)+y_b_u(z1))/d1];
                        v_approach_2=[vx_drop(j,i) vy_drop(j,i)]*[(-x_drop(j,i)+x_b(z2))/d2; (-y_drop(j,i)+y_b_l(z2))/d2];
                        
                        c1=x_b(z1);
                        c2=x_b(z2);
                        c3=y_b_u(z1);
                        c4=y_b_l(z2);
                        
                        if d1<=scale_b_2 && v_approach_1>0
                            if d1<=scale_b_1
                                %                         X component equation
                                D(jc,jc)=D(jc,jc)+((c1-x_drop(j,i))/d1)*(kb/(d1-R))*(((x_drop(j,i)-c1))/d1);
                                D(jc,counter+jc)=D(jc,counter+jc)+((c3-y_drop(j,i))/d1)*(kb/(d1-R))*(((x_drop(j,i)-c1))/d1);
                                %                         Y component equation
                                D(counter+jc,jc)=D(counter+jc,jc)+((c1-x_drop(j,i))/d1)*(kb/(d1-R))*(((y_drop(j,i)-c3))/d1);
                                D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+((c3-y_drop(j,i))/d1)*(kb/(d1-R))*(((y_drop(j,i)-c3))/d1);
                            else
                                D(counter+jc,jc)=D(counter+jc,jc)+((c1-x_drop(j,i))/d1)*(kb/(d1-R))*(((y_drop(j,i)-c3))/d1);
                                D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+((c3-y_drop(j,i))/d1)*(kb/(d1-R))*(((y_drop(j,i)-c3))/d1);
                            end
                        end
                        
                        if d2<=scale_b_2 && v_approach_2>0
                            if d2<=scale_b_1
                                %                         X component equation
                                D(jc,jc)=D(jc,jc)+((c2-x_drop(j,i))/d2)*(kb/(d2-R))*(((x_drop(j,i)-c2))/d2);
                                D(jc,counter+jc)=D(jc,counter+jc)+((c4-y_drop(j,i))/d2)*(kb/(d2-R))*(((x_drop(j,i)-c2))/d2);
                                %                         Y component equation
                                D(counter+jc,jc)=D(counter+jc,jc)+((c2-x_drop(j,i))/d2)*(kb/(d2-R))*(((y_drop(j,i)-c4))/d2);
                                D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+((c4-y_drop(j,i))/d2)*(kb/(d2-R))*(((y_drop(j,i)-c4))/d2);
                            else
                                D(counter+jc,jc)=D(counter+jc,jc)+((c2-x_drop(j,i))/d2)*(kb/(d2-R))*(((y_drop(j,i)-c4))/d2);
                                D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+((c4-y_drop(j,i))/d2)*(kb/(d2-R))*(((y_drop(j,i)-c4))/d2);
                            end
                        end
                    else
                        
                        % interaction for the drops in the straight channel
                        d1=y_ec_u-y_drop(j,i);
                        d2=y_drop(j,i)-y_ec_l;
                        v_approach=vy_drop(j,i);
                        
                        c3=y_ec_u;
                        c4=y_ec_l;
                        
                        if d1<=scale_b_1 && v_approach>0
                            D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+(kd/(d1-R))*(((y_drop(j,i)-c3))/d1);
                        end
                        
                        if d2<=scale_b_1 && v_approach<0
                            D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+(kd/(d2-R))*(((y_drop(j,i)-c4))/d2);
                        end
                    end
                end
                
                
                if c_old_1-c+p2+1 ~= c_old_1
                    for ci = 1:(c-p2-1)
                        d_sol_big(:,c_old_1-c+p2+1+ci)=zeros(size(d_sol_big,1),1);
                    end
                end
                
                
                V_final = [v_sol_big g_sol_big];
                D_final = [d_sol_big D];
                
                A = [V_final; D_final];
                B= [c_sol_big; f];
                
                % Solving for the new velocities
                
                % Re-scaling the matrices to avoid instability
                scale_A = zeros(size(A,1),1);
                A_scaled = zeros(size(A));
                B_scaled = zeros(size(B));
                for var1 = 1:size(A,1)
                    sum_scale=0;
                    for var2 = 1:size(A,2)
                        sum_scale = sum_scale + abs(A(var1,var2));
                    end
                    power_2 = double(-int16(log2(sum_scale))-1);
                    scale_A(var1,1) =  2^power_2;
                    A_scaled(var1,:)= A(var1,:).*scale_A(var1,1);
                    B_scaled(var1,:) = B(var1,:).*scale_A(var1,1);
                end
                
                U_scaled = A_scaled\B_scaled;
                U = U_scaled;
                
                vel_start_index = size(c_sol_big,1)+1;
                guess(1:counter)=U(vel_start_index:vel_start_index+counter-1,1);
                guess(counter+1:2*counter)=U(vel_start_index+counter:vel_start_index+2*counter-1,1);
                
                vx_drop(drops_test,i)=guess(1:counter);
                vy_drop(drops_test,i)=guess(counter+1:2*counter);
                
                v1_d=0;
                for jc_d=1:counter
                    
                    j_d=drops_test(jc_d);
                    c_d=2*no(jc_d)+2;
                    
                    yfd = zeros(c_d,2);
                    yfd(:,1) = y_f_d_big(v1_d+1:v1_d+c_d,1);
                    yfd(1,2) =0;
                    yfd(2:c_d-1,2) = [sqrt(vx_drop(ino(1:no(jc_d),jc_d,1),i).^2+vy_drop(ino(1:no(jc_d),jc_d,1),i).^2);...
                        sqrt(vx_drop(ino(1:no(jc_d),jc_d,1),i).^2+vy_drop(ino(1:no(jc_d),jc_d,1),i).^2)];
                    yfd(1,c_d)=0;
                    
                    b=0;
                    for kk=2:2:c_d-2
                        b=b+(yfd(kk,1)-yfd(kk+1,1))*yfd(kk,2);
                    end
                    
                    csol_d = [yfd(:,2); flowrate-b];
                    vsol_d = v_sol_big(v1_d+1:v1_d+c_d+1,v1_d+1:v1_d+c_d+1);
                    coeff_sol = vsol_d\csol_d;
                    
                    
                    
                    p1_d=find(yfd(:,1)==y_drop(j_d,i)+R);
                    p2_d=find(yfd(:,1)==y_drop(j_d,i)-R);
                    
                    p1_d=p1_d(1);
                    p2_d=p2_d(1);
                    
                    A=coeff_sol(1);
                    B1=coeff_sol(p1_d);
                    C1=coeff_sol(p1_d+1);
                    B2=coeff_sol(p2_d+1);
                    C2=coeff_sol(p2_d+2);
                    
                    v1=(A/3*((yfd(p1_d-1,1)^3)-(yfd(p1_d,1)^3))+B1/2*((yfd(p1_d-1,1)^2)-(yfd(p1_d,1)^2))+C1*(yfd(p1_d-1,1)-yfd(p1_d,1)))/(yfd(p1_d-1,1)-yfd(p1_d,1));
                    v2=(A/3*((yfd(p2_d,1)^3)-(yfd(p2_d+1,1)^3))+B2/2*((yfd(p2_d,1)^2)-(yfd(p2_d+1,1)^2))+C2*(yfd(p2_d,1)-yfd(p2_d+1,1)))/(yfd(p2_d,1)-yfd(p2_d+1,1));
                    
                    % Making of the f matrix [DU=f]
                    f_d(jc_d)=-kf*(v1+v2)/2;
                    f_d(counter+jc_d)=-kf*k1*(v1-v2);
                    
                    v1_d = v1_d+ c_d+1;
                    
                end
                
                % Solving for the new velocities
                U_d = D\f_d ;
                
                
                % Finding the error in the solutions
                err=U_d-[vx_drop(drops_test,i);vy_drop(drops_test,i)];

                % Updating the guess
                guess=[vx_drop(drops_test,i);vy_drop(drops_test,i)]+alpha_value*err;
                
                % relative error as a test of convergence
                rel_err=norm(err/vx_guess);
                rer(ii,i)=rel_err;
                
                if ii>1
                    
                    err_rate=abs(rer(ii)-rer(ii-1));
                end
                
                ii=ii+1;
                % Escape pod in terms of no convergence
                lim(i,p)=lim(i,p)+1;
                if lim(i,p)==1e3
                    rel_err=0;
                    disp('Error in time step:')
                    disp(i)
                end
                
            end
            
            vx_drop(drops_test,i)=U_d(1:counter);
            vy_drop(drops_test,i)=U_d(counter+1:2*counter);
            
        end
        
        %POSITION UPDATE
        x_drop(:,i+1)=x_drop(:,i)+vx_drop(:,i)*t_step;
        y_drop(:,i+1)=y_drop(:,i)+vy_drop(:,i)*t_step;
        
        % GUESS UPDATE
        v_guess=[vx_drop(:,i); vy_drop(:,i)];
        
        
        if rem(i,100)==0
            disp('i')
            disp(i)
        end
        sim_time(i)=toc;
    end
    
    ts=find(drop_in_test==max(drop_in_test), 1 );
    tf=t;
    
    t_sim_steps(p)=tf-ts;
    t_sim(p)=sim_time(tf)-sim_time(ts);
    
    xd(:,:,p)=x_drop;
    yd(:,:,p)=y_drop;
    vxd(:,:,p)=vx_drop;
    vyd(:,:,p)=vy_drop;
    
    %     save(name)
end

%%
ima=round(linspace(1,t,5));
d=ceil((length(ima))/2);

[x,y]=cylinder(R,100);
cc=zeros(n,3);
a=[1 0 0];
b=[0 0.5 0.5];

for ii=1:n
    cc(ii,:)=b;
end

for p=1:length(L0)
    for ii=1:length(ima)
        figure;
        q=ima(ii);
        
        xch1=[-4*R x_b max(x_b)+4*R max(x_b)+4*R -4*R];
        ych1=[y_ec_u y_b_u y_ec_u y_ec_u+max(y_b_u) y_ec_u+max(y_b_u)];
        xch2=xch1;
        ych2=[y_ec_l y_b_l y_ec_l y_ec_l+min(y_b_l) y_ec_l+min(y_b_l)];
        
        fill(xch1,ych1,'k')
        hold on
        fill(xch2,ych2,'k')
        hold on
        
        for k=drop_in_test_initial(q,p):drop_in_test_final(q,p)
            
            fill(x(1,:)+xd(k,q,p),y(1,:)+yd(k,q,p),cc(k,:))
            hold on
            axis equal
        end
        hold off
        axis equal
        set(gca, 'XTick', [],'YTick', []);
        axis([min(x_b)-3*R max(x_b)+0*R min(y_b_l)-0.4*R max(y_b_u)+0.4*R])
    end
end