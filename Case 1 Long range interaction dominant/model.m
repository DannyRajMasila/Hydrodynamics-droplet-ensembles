function force = model(r)

global Y
global n
global K K0

id=1.195;

x_drop=r(1:n,1);
y_drop=r(n+1:2*n,1);

x_ij=zeros(n);
y_ij=zeros(n);
r_ij=zeros(n);
e_ijx=zeros(n);
e_ijy=zeros(n);

for p=1:n
    for q=1:n
        x_ij(p,q)=x_drop(p)-x_drop(q);
        y_ij(p,q)=y_drop(p)-y_drop(q);
        r_ij(p,q)=norm([x_ij(p,q), y_ij(p,q)]);
        e_ijx(p,q)=x_ij(p,q)/r_ij(p,q);
        e_ijy(p,q)=y_ij(p,q)/r_ij(p,q);
    end
end

%     r_ij;
drift_term_x = zeros(n,1);
dipolar_term_x = zeros(n,1);
adhesive_term_x = zeros(n,1);
repulsive_term_x = zeros(n,1);

drift_term_y = zeros(n,1);
dipolar_term_y = zeros(n,1);
adhesive_term_y = zeros(n,1);
repulsive_term_y = zeros(n,1);

for t=1:n
    for s=1:n
        if t~=s
            drift_term_x(t,1)=drift_term_x(t,1)+ (1/r_ij(t,s))^2;
            dipolar_term_x(t,1)=dipolar_term_x(t,1)+(-2*((1/r_ij(t,s))^2)*(e_ijx(t,s)*e_ijx(t,s)));
            adhesive_term_x(t,1)=adhesive_term_x(t,1)-Y*((1/(r_ij(t,s)-id))^2)*e_ijx(t,s);
            repulsive_term_x(t,1)=repulsive_term_x(t,1)+0.5*K*((1/(r_ij(t,s)-id))^12)*e_ijx(t,s);
            
            drift_term_y(t,1)=0;
            dipolar_term_y(t,1)=dipolar_term_y(t,1)+(-2*((1/r_ij(t,s))^2)*(e_ijy(t,s)*e_ijx(t,s)));
            adhesive_term_y(t,1)=adhesive_term_y(t,1)-Y*((1/(r_ij(t,s)-id))^2)*e_ijy(t,s);
            repulsive_term_y(t,1)=repulsive_term_y(t,1)+0.5*K*((1/(r_ij(t,s)-id))^12)*e_ijy(t,s);
        end
    end
end
x = K0*(drift_term_x+dipolar_term_x)+adhesive_term_x+repulsive_term_x;
y = K0*(drift_term_y+dipolar_term_y)+adhesive_term_y+repulsive_term_y;

force=[x,y];
end