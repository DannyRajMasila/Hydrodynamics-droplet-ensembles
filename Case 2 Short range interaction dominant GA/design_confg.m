
function obj_fn = design_confg(x_drop,y_drop,start_var,end_var,~,ObjFn)
global R

start_d=start_var;
end_d=end_var;

Z=[x_drop(start_d:end_d),y_drop(start_d:end_d)];
cx = mean(Z(:,1));
cy = mean(Z(:,2));

D = sqrt((Z(:,1)-cx).^2+(Z(:,2)-cy).^2);
Y = sort(D);

obj_fn= (1/R^2).*sum(var(Y-ObjFn));
end



