function [out] = FVyz(T,dY,dZ)

ou  = FVyGrad(FVzGrad(T,dZ),dY);
out = ou(2:end-1,2:end-1);

end