

function var_c = MAS_solver_AD_Optimization(start_var,end_var,n,variables,code,ObjFn,~)

parameters_final

global flowrate
global x_b
global y_b_u
global y_b_l
global y_ec_u
global y_ec_l
global alpha_value
global lambda_value
global gamma_value
global R


L0=1;
R1=1.25e-4; % Half of inlet channel dimension
x0=0; %Initial position of drop1 wrt 'x'
y0=0; %Initial position of drop1 wrt 'y'
t=100000;
t_step= 1e-4; % single time step
% Need for a smaller time step than Danny's code

vx_guess=0.11;
vy_guess=0.01;
y_ec_u=1*R1;
y_ec_l=-1*R1;

L = zeros(1,n);
for var_y=1:n-1
    L(var_y)=variables(var_y);
end

Q1= 400;
Q = 40; % Fix Q

alp = variables(var_y+1);
H = variables(var_y+2);
%--------------------------------------------------------------------

drop_in_test_initial=zeros(t,length(L0));
drop_in_test_final=zeros(t,length(L0));
vxd=zeros(n,t,length(L0));
vyd=zeros(n,t,length(L0));
xd=zeros(n,t+1,length(L0));
yd=zeros(n,t+1,length(L0));

% % NEW GEOMETRY
x_b=linspace(0,W,1002);
x_b_1 = linspace(0,alp*W,ceil(alp*1002));
y_b_u_1=(heaviside(x_b_1+R1)-heaviside(x_b_1-alp*W)).*((h*alp*W-h*x_b_1+H*x_b_1)/(alp*W))...
   +(heaviside(x_b_1-alp*W)-heaviside(x_b_1-1.2*W)).*((H*W-h*alp*W+h*x_b_1-H*x_b_1)/(W-alp*W));

x_b_2 = linspace(alp*W,W,1002-ceil(alp*(1002)));
slope_small = 0.01;
y_b_u_2= slope_small.*x_b_2+(y_b_u_1(:,end)-slope_small*alp*W);

y_b_u = [y_b_u_1 y_b_u_2];
y_b_l=-y_b_u;

% ADDING ROUGHNESS TO WALL
% Asymmetry in the top-down roughness
mul1= 1; mul2= 0.1;
y_b_u=y_b_u+0.001*mul1*R1*rand(1,length(y_b_u));
y_b_l=y_b_l-0.001*mul2*R1*rand(1,length(y_b_l));

y_b_u(1)=y_ec_u;
y_b_l(1)=y_ec_l;
y_b_u(length(y_b_u))=y_ec_u;
y_b_l(length(y_b_l))=y_ec_l;

for p=1:length(L0)
    
    Q2=Q;
    flowrate=((Q1+Q2)*1e-6*1e-3)./(60*(2*R1));
    
    %initialization
    vx_drop=zeros(n,t);
    vy_drop=zeros(n,t);
%     vd_drop=zeros(n,t);
    x_drop=zeros(n,t+1);
    y_drop=zeros(n,t+1);
    
    % NEW SPACING
    
    x_drop(1,1) = x0;
    y_drop(1,1) = y0;
    for k=2:n
        x_drop(k,1)=x_drop(k-1,1)-L(1,k-1)*R;
        y_drop(k,1)=y0;
    end
    
    %INITIAL VELOCITIES-GUESS
    v_guess=[vx_guess*ones(n,1);vy_guess*ones(n,1)];
    
    %VELOCITY IN THE INLET AND OUTLET CHANNELS
    v_inletoutlet=1e-9*(Q1+Q2)/((y_ec_u-y_ec_l)*(2*R1))/60;
    
    % %QUASI STEADY STATE SIMULATION------------------------------------------
    drop_in_test=zeros(t,1);
    rer=zeros(100,t);
    
    for i=1:t
        
        if (x_drop(end,i)<alp*W+3*R)
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
            
            % DROP INDEX
            drop_index = zeros(drops_test(end),1);
            drop_index(drops_test,1)= (1:size(drops_test,1))';
            
            
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
                
                % parameters
                
                no=zeros(counter,1);
                ino=zeros(10,counter,3);
                inb=zeros(counter,2);
                
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

                while (rel_err>1e-8) && (abs(err_rate)>0.2*1e-10)
                    
                    % Linear code solved iteratively along with the non-linear
                    % code
                    alpha_value=1; %#ok<NASGU>
                    lambda_value=0; %#ok<NASGU>
                    gamma_value=0; %#ok<NASGU>
                    [~,~,guess_lin] = linear_code(drop_index,guess,counter,drops_test,no,ino,inb,x_drop(:,i),y_drop(:,i),vx_drop(:,i),vy_drop(:,i));
                    
                    alpha_value=0;
                    lambda_value=1;
                    gamma_value=1;
                    
                    [U_d,err,guess] = non_linear_code(guess_lin,guess, counter,drops_test,no,ino,inb,x_drop(:,i),y_drop(:,i),vx_drop(:,i),vy_drop(:,i));
                    
                    % relative error as a test of convergence
                    rel_err=norm(err/vx_guess);
                    rer(ii,i)=rel_err;
                    
                    if ii>1
                        
                        err_rate=rer(ii)-rer(ii-1);
                    end
                    
                    ii=ii+1;
                end
                
                vx_drop(drops_test,i)=U_d(1:counter);
                vy_drop(drops_test,i)=U_d(counter+1:2*counter);
                
            end
            
            %POSITION UPDATE
            x_drop(:,i+1)=x_drop(:,i)+vx_drop(:,i)*t_step;
            y_drop(:,i+1)=y_drop(:,i)+vy_drop(:,i)*t_step;
            
            % GUESS UPDATE
            v_guess=[vx_drop(:,i); vy_drop(:,i)];
            
        else
            break;
        end
    end
    
    xd(:,:,p)=x_drop;
    yd(:,:,p)=y_drop;
    vxd(:,:,p)=vx_drop;
    vyd(:,:,p)=vy_drop;
    t=i-1;
end

var_c = feval(@design_confg,x_drop(:,t),y_drop(:,t),start_var,end_var,code,ObjFn);
% disp([variables,var_c])
end

