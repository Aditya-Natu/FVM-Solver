%Code developed by Aditya Natu%

function du_dz = d_dz(U,dZ)

U_N = U(2:end-2,:);
U_S = U(3:end-1,:);

dZ_ext = [dZ(:,1),dZ,dZ(:,end)];

dZ_N = dZ_ext(1:end-1,:);
dZ_S = dZ_ext(2:end-0,:);

U_mid = (U_N.*dZ_S + U_S.*dZ_N)./(dZ_S + dZ_N);

U_No = [U(1,:  );U_mid];
U_So = [U_mid;U(end,:)];

du_dz_mid = (U_So - U_No)./dZ_ext;

du_dz_s = 2*(U(end,:) - U(end-1,:))./(dZ_ext(end,:));
du_dz_n = 2*(U(2  ,:) - U(1    ,:))./(dZ_ext(1  ,:));

du_dz = [du_dz_n;du_dz_mid;du_dz_s];

end
