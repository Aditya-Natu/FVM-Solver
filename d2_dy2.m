%Developed and tested by Aditya Natu%

function [sol] = d2_dy2(T,dY)

U_w = T(:,1:end-2);
U_e = T(:,3:end  );
U_p = T(:,2:end-1);

dY_ext = [dY(1,:);dY;dY(end,:)];

dY_W = [0*dY_ext(:,1  ),dY_ext(:,1:end-1)];
dY_E = [dY_ext(:,2:end-0),0*dY_ext(:,end)];

WF = 2*(U_p - U_w)./(dY_W + dY_ext);
EF = 2*(U_e - U_p)./(dY_E + dY_ext);

d2U_dY2 = (EF - WF)./dY_ext;

dU_dY_Wp = ((T(:,2    ).*dY_ext(:,2    ) + T(:,3    ).*dY_ext(:,1  ))./(dY_ext(:,1    )+dY_ext(:,2  )) - T(:,1  ))./dY_ext(:,1  );
d2U_dY2_WF =  2*(dU_dY_Wp - WF(:,1  ))./dY_ext(:,1  );

dU_dY_Ep = (T(:,end) - (T(:,end-1).*dY_ext(:,end-1) + T(:,end-2).*dY_ext(:,end))./(dY_ext(:,end-1)+dY_ext(:,end)))./dY_ext(:,end);
d2U_dY2_EF =  2*(EF(:,end) - dU_dY_Ep)./dY_ext(:,end);

sol = [d2U_dY2_WF, d2U_dY2, d2U_dY2_EF];

end
