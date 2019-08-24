function dr = position(~,r)
global n
dr = zeros(2*n,1);

F= model(r);

for i=1:2*n
dr(i)=F(i);
end

end