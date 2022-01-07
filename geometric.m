%% codegen
function sol = geometric(a,b,n,r)
t1 = (b-a).*(1-r)./(1 - r.^(n-1));
m = 0:1:n-1;
sol = a + t1.*(1-r.^m)./(1-r);
end