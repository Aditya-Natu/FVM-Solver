%Code developed by Aditya Natu%

function du_dy = d_dy(U,dY)

U_W = U(:,2:end-2);
U_E = U(:,3:end-1);

dY_ext = [dY(1,:);dY;dY(end,:)];

dY_W = dY_ext(:,1:end-1);
dY_E = dY_ext(:,2:end-0);

U_mid = (U_W.*dY_E + U_E.*dY_W)./(dY_E + dY_W);

U_We = [U(:,1  ),U_mid];
U_Ea = [U_mid,U(:,end)];

du_dy_mid = (U_Ea - U_We)./dY_ext;

du_dy_e = 2*(U(:,end) - U(:,end-1))./(dY_ext(:,end));
du_dy_w = 2*(U(:,2  ) - U(:,1    ))./(dY_ext(:,1  ));

du_dy = [du_dy_w,du_dy_mid,du_dy_e];

end
