function [div_Y] = FVy2Grad(T,dY)

U_w = T(2:end-1,1:end-2);
U_e = T(2:end-1,3:end-0);
U_p = T(2:end-1,2:end-1);

dY_W = [0*dY(:,1  ),dY(:,1:end-1)];
dY_E = [dY(:,2:end-0),0*dY(:,end)];

WF = 2*(U_w - U_p)./(dY_W + dY);
EF = 2*(U_e - U_p)./(dY_E + dY);

div_Y = (EF + WF)./dY;

end