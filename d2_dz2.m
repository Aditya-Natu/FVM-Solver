%Developed and tested by Aditya Natu%

function [sol] = d2_dz2(T,dZ)

U_n = T(1:end-2,:);
U_s = T(3:end  ,:);
U_p = T(2:end-1,:);

dZ_ext = [dZ(:,1),dZ,dZ(:,end)];

dZ_N = [0*dZ_ext(1  ,:);dZ_ext(1:end-1,:)];
dZ_S = [dZ_ext(2:end-0,:);0*dZ_ext(end,:)];

NF = 2*(U_p - U_n)./(dZ_N + dZ_ext);
SF = 2*(U_s - U_p)./(dZ_S + dZ_ext);

d2U_dZ2 = (SF - NF)./dZ_ext;

dU_dZ_Np = ((T(2    ,:).*dZ_ext(2    ,:) + T(3    ,:).*dZ_ext(1  ,:))./(dZ_ext(1    ,:)+dZ_ext(2  ,:)) - T(1  ,:))./dZ_ext(1  ,:);
d2U_dZ2_NF =  2*(dU_dZ_Np - NF(1  ,:))./dZ_ext(1  ,:);

dU_dZ_Sp = (T(end,:) - (T(end-1,:).*dZ_ext(end-1,:) + T(end-2,:).*dZ_ext(end,:))./(dZ_ext(end-1,:)+dZ_ext(end,:)))./dZ_ext(end,:);
d2U_dZ2_SF =  2*(SF(end,:) - dU_dZ_Sp)./dZ_ext(end,:);

sol = [d2U_dZ2_NF; d2U_dZ2; d2U_dZ2_SF];

end
