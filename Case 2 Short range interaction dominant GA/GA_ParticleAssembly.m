% GA code for the identifying the conditions to make desired assemblies of
% droplets

clc
clear

% Shapes for which desired features has been given: 
% 'circle5', 'circle6', 'triangle6', 'square8'

shape='circle6';
OBJ_fun

tstart=tic;
parameters_final

% Number of preambles required: Preamble droplets are inert droplets that
% participate in the pattern formation and not the fusion into particle.
% Also, in the current formulation, preamble droplets always preceed the
% active ones.
no_preamble=1;

% Total number of droplets in the simulation
n=size(ObjFn,1)+no_preamble;

% Bounds
LB = [3.*ones(n-1,1);0.3;0.0025];
UB = [10.*ones(n-1,1);0.3;0.0025];

start_var =no_preamble+1;
end_var = n;
code =13;
nvar = (n-1)+2;

% GA solver inbuilt- run in multiple parallel cores
options = gaoptimset('UseParallel','always','Display','iter');
[X, fval] = ga(@(variables) MAS_solver_AD_Optimization(start_var,end_var,n,variables,code,ObjFn,tstart),nvar,[],[],[],[],LB,UB,[],options);

disp(X)