% System parameters
global R
mu=0.019; %viscosity cp
R=1.09e-4; %radius of the drop in m 


% Tuned parameters
kf=6*pi*mu*R;
k1=0.25;
beta_value=1;
kd=3*6*pi*mu*R^2;
kb=5*6*pi*mu*R^2;
scale_d=2.4*R;
scale_b_1=1.01*R;
scale_b_2=12*R;
R1=1.25e-4; % Half of inlet channel dimension


h=R1;
W=21*2*R1;
alp = 0.3;
H = 0.0025;

