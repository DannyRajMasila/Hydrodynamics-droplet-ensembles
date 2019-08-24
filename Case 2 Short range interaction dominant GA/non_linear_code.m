function [U,err1,guess_new] = non_linear_code(guess_lin,guess,counter,drops_test,no,ino,inb, x_drop,y_drop,vx_drop,vy_drop)
global flowrate
global x_b
global y_b_u
global y_b_l
global y_ec_u
global y_ec_l
global alpha_value
global lambda_value
global gamma_value

parameters_final
vx_drop(drops_test)=guess_lin(1:counter);
vy_drop(drops_test)=guess_lin(counter+1:2*counter);

D=zeros(2*counter,2*counter);
f=zeros(2*counter,1);

for jc=1:counter
    j=drops_test(jc);
    %     m=x_drop(j);
    c=2*no(jc)+2;
    
    y_f=zeros(c,2);
    
    y_f(1,1)=inb(jc,1);
    y_f(1,2)=0;
    y_f(c,1)=inb(jc,2);
    y_f(c,2)=0;
    y_f(2:c-1,1)=[ino(1:no(jc),jc,2) ; ino(1:no(jc),jc,3)];
    y_f(2:c-1,2)=[sqrt(vx_drop(ino(1:no(jc),jc,1)).^2+vy_drop(ino(1:no(jc),jc,1)).^2);...
        sqrt(vx_drop(ino(1:no(jc),jc,1)).^2+vy_drop(ino(1:no(jc),jc,1)).^2)];
    
    [a, ix]=sort(y_f(:,1),'descend');
    temp=zeros(c,1);
    for q=1:c
        temp(q)=y_f(ix(q),2);
    end
    y_f(:,1)=a;
    y_f(:,2)=temp;
    
    %Constructing the Flow solution Matrix (v_sol); Size=(c+1)x(c+1)
    v_sol=zeros(c+1,c+1);
    
    %'y and 1' terms in c of c+1 rows
    for g=1:c
        hh=2;
        gg=ceil(g/2)*hh;
        for h=1:2
            v_sol(g,gg)=y_f(g,1)^(2-h);
            gg=gg+1;
        end
    end
    
    %First column
    v_sol(1,1)=y_f(1,1)^2;
    v_sol(c,1)=y_f(c,1)^2;
    v_sol(1:c,1)=y_f(:,1).^2;
    
    %last row
    for h=2:2:c
        v_sol(c+1,h)=(y_f(h-1,1)^2-y_f(h,1)^2)/2;
    end
    for h=3:2:c+1
        v_sol(c+1,h)=y_f(h-2,1)-y_f(h-1,1);
    end
    for h=1:2:c
        v_sol(c+1,1)=v_sol(c+1,1)+(y_f(h,1)^3-y_f(h+1,1)^3)/3;
    end
    b=0;
    for kk=2:2:c-2
        b=b+(y_f(kk,1)-y_f(kk+1,1))*y_f(kk,2);
    end
    c_sol=[y_f(:,2);(flowrate)-b];
    
    %---------------------------------------------------------------------
    % Solution for the flow
    coeff_sol=v_sol\c_sol;
    %---------------------------------------------------------------------
    
    %AVERAGE VELOCITIES
    p1=find(y_f(:,1)==y_drop(j)+R);
    p2=find(y_f(:,1)==y_drop(j)-R);
    
    p1=p1(1);
    p2=p2(1);
    
    A=coeff_sol(1);
    B1=coeff_sol(p1);
    C1=coeff_sol(p1+1);
    B2=coeff_sol(p2+1);
    C2=coeff_sol(p2+2);
    
    v1=(A/3*((y_f(p1-1,1)^3)-(y_f(p1,1)^3))+B1/2*((y_f(p1-1,1)^2)-(y_f(p1,1)^2))+C1*(y_f(p1-1,1)-y_f(p1,1)))/(y_f(p1-1,1)-y_f(p1,1));
    v2=(A/3*((y_f(p2,1)^3)-(y_f(p2+1,1)^3))+B2/2*((y_f(p2,1)^2)-(y_f(p2+1,1)^2))+C2*(y_f(p2,1)-y_f(p2+1,1)))/(y_f(p2,1)-y_f(p2+1,1));
    
    % Making of the f matrix [DU=f]
    f(jc)=-kf*(v1+v2)/2;
    f(counter+jc)=-kf*k1*(v1-v2);
    
    % Assembling the D matrix [DU=f]
    D(jc,jc)=-beta_value*kf;
    D(counter+jc,counter+jc)=-beta_value*kf;
    
    % Effect of drop-drop interactions incorporated in D [DU=f]
    
    for jl=1:counter
        l=drops_test(jl);
        d=sqrt((x_drop(j)-x_drop(l))^2+(y_drop(j)-y_drop(l))^2);
        
        v_approach=([vx_drop(j) vy_drop(j)]*[x_drop(l)-x_drop(j);  y_drop(l)-y_drop(j)]...
            +[vx_drop(l) vy_drop(l)]*[-x_drop(l)+x_drop(j) ; -y_drop(l)+y_drop(j)]);
        if d==0 || d>=scale_d || v_approach<0    % APPROACH VELOCITY POSITIVE FOR INTERACTION TO OCCUR- NEEDS VERIFICATION
        else
            %                     X component equations
            D(jc,jc)=D(jc,jc)+((x_drop(l)-x_drop(j))/d)*(kd/(d-2*R))*(((x_drop(j)-x_drop(l)))/d);
            D(jc,jl)=D(jc,jl)+((x_drop(j)-x_drop(l))/d)*(kd/(d-2*R))*(((x_drop(j)-x_drop(l)))/d);
            D(jc,counter+jc)=D(jc,counter+jc)+((y_drop(l)-y_drop(j))/d)*(kd/(d-2*R))*(((x_drop(j)-x_drop(l)))/d);
            D(jc,counter+jl)=D(jc,counter+jl)+((y_drop(j)-y_drop(l))/d)*(kd/(d-2*R))*(((x_drop(j)-x_drop(l)))/d);
            %                     Y component equations
            D(counter+jc,jc)=D(counter+jc,jc)+((x_drop(l)-x_drop(j))/d)*(kd/(d-2*R))*(((y_drop(j)-y_drop(l)))/d);
            D(counter+jc,jl)=D(counter+jc,jl)+((x_drop(j)-x_drop(l))/d)*(kd/(d-2*R))*(((y_drop(j)-y_drop(l)))/d);
            D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+((y_drop(l)-y_drop(j))/d)*(kd/(d-2*R))*(((y_drop(j)-y_drop(l)))/d);
            D(counter+jc,counter+jl)=D(counter+jc,counter+jl)+((y_drop(j)-y_drop(l))/d)*(kd/(d-2*R))*(((y_drop(j)-y_drop(l)))/d);
        end
    end
    
    % Effect of drop-boundary interactions incorporated in D [DU=f]
    
    if x_drop(j)>=min(x_b) && x_drop(j)<=max(x_b)
        % interaction for the drops in the 2D section
        
        [d1, z1]=min((x_b-x_drop(j)).^2+(y_b_u-y_drop(j)).^2);
        [d2, z2]=min((x_b-x_drop(j)).^2+(y_b_l-y_drop(j)).^2);
        d1=sqrt(d1);
        d2=sqrt(d2);
        v_approach_1=[vx_drop(j) vy_drop(j)]*[(-x_drop(j)+x_b(z1))/d1; (-y_drop(j)+y_b_u(z1))/d1];
        v_approach_2=[vx_drop(j) vy_drop(j)]*[(-x_drop(j)+x_b(z2))/d2; (-y_drop(j)+y_b_l(z2))/d2];
        
        c1=x_b(z1);
        c2=x_b(z2);
        c3=y_b_u(z1);
        c4=y_b_l(z2);
        
        if d1<=scale_b_2 && v_approach_1>0
            if d1<=scale_b_1
                %                         X component equation
                D(jc,jc)=D(jc,jc)+((c1-x_drop(j))/d1)*(kb/(d1-R))*(((x_drop(j)-c1))/d1);
                D(jc,counter+jc)=D(jc,counter+jc)+((c3-y_drop(j))/d1)*(kb/(d1-R))*(((x_drop(j)-c1))/d1);
                %                         Y component equation
                D(counter+jc,jc)=D(counter+jc,jc)+((c1-x_drop(j))/d1)*(kb/(d1-R))*(((y_drop(j)-c3))/d1);
                D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+((c3-y_drop(j))/d1)*(kb/(d1-R))*(((y_drop(j)-c3))/d1);
            else
                D(counter+jc,jc)=D(counter+jc,jc)+((c1-x_drop(j))/d1)*(kb/(d1-R))*(((y_drop(j)-c3))/d1);
                D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+((c3-y_drop(j))/d1)*(kb/(d1-R))*(((y_drop(j)-c3))/d1);
            end
        end
        
        if d2<=scale_b_2 && v_approach_2>0
            if d2<=scale_b_1
                %                         X component equation
                D(jc,jc)=D(jc,jc)+((c2-x_drop(j))/d2)*(kb/(d2-R))*(((x_drop(j)-c2))/d2);
                D(jc,counter+jc)=D(jc,counter+jc)+((c4-y_drop(j))/d2)*(kb/(d2-R))*(((x_drop(j)-c2))/d2);
                %                         Y component equation
                D(counter+jc,jc)=D(counter+jc,jc)+((c2-x_drop(j))/d2)*(kb/(d2-R))*(((y_drop(j)-c4))/d2);
                D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+((c4-y_drop(j))/d2)*(kb/(d2-R))*(((y_drop(j)-c4))/d2);
            else
                D(counter+jc,jc)=D(counter+jc,jc)+((c2-x_drop(j))/d2)*(kb/(d2-R))*(((y_drop(j)-c4))/d2);
                D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+((c4-y_drop(j))/d2)*(kb/(d2-R))*(((y_drop(j)-c4))/d2);
            end
        end
    else
        
        % interaction for the drops in the straight channel
        d1=y_ec_u-y_drop(j);
        d2=y_drop(j)-y_ec_l;
        v_approach=vy_drop(j);
        
        c3=y_ec_u;
        c4=y_ec_l;
        
        if d1<=scale_b_1 && v_approach>0
            D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+(kd/(d1-R))*(((y_drop(j)-c3))/d1);
        end
        
        if d2<=scale_b_1 && v_approach<0
            D(counter+jc,counter+jc)=D(counter+jc,counter+jc)+(kd/(d2-R))*(((y_drop(j)-c4))/d2);
        end
    end
end

% Solving for the new velocities
% Re-scaling the matrices to avoid instability
scale_D = zeros(size(D,1),1);
D_scaled = zeros(size(D));
f_scaled = zeros(size(f));
for var1 = 1:size(D,1)
    sum_scale=0;
    for var2 = 1:size(D,2)
        sum_scale = sum_scale + abs(D(var1,var2));
    end
    power_2 = double(-int16(log2(sum_scale))-1);
    scale_D(var1,1) =  2^power_2;
    D_scaled(var1,:)= D(var1,:).*scale_D(var1,1);
    f_scaled(var1,:) = f(var1,:).*scale_D(var1,1);
end

U_scaled = D_scaled\f_scaled;
U = U_scaled;

% Finding the error in the solutions
err1=U-[guess(1:counter);guess(counter+1:2*counter)];
err2=U-[vx_drop(drops_test);vy_drop(drops_test)];
err3=[vx_drop(drops_test);vy_drop(drops_test)]-[guess(1:counter);guess(counter+1:2*counter)];


% Updating the guess
guess_new=[guess(1:counter);guess(counter+1:2*counter)]+alpha_value*err1+...
                                        gamma_value*err3+lambda_value*err2;


end
