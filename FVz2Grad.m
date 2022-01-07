function [div_Z] = FVz2Grad(T,dZ)

V_n = T(1:end-2,2:end-1);
V_s = T(3:end-0,2:end-1);
V_p = T(2:end-1,2:end-1);

dZ_N = [0*dZ(1  ,:);dZ(1:end-1,:)];
dZ_S = [dZ(2:end-0,:);0*dZ(end,:)];

NF = 2*(V_n - V_p)./(dZ_N + dZ);
SF = 2*(V_s - V_p)./(dZ_S + dZ);

div_Z = (SF + NF)./dZ;

end