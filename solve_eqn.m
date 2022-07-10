function [U_ext] = solve_eqn(Cy,Cz,Dy,Dz,Prp,RHS,dY,dZ,NB,EB,WB,SB)

[sz,sy] = size(dY);

dYE = [dY(:,2:end),  zeros(sz,1)];
dYW = [zeros(sz,1),dY(:,1:end-1)];
dZN = [zeros(1,sy);dZ(1:end-1,:)];
dZS = [dZ(2:end,:);  zeros(1,sy)];

CyE = (Cy(2:end-1,2:end-1).*dYE + Cy(2:end-1,3:end-0).*dY)./(dY + dYE);
CyW = (Cy(2:end-1,2:end-1).*dYW + Cy(2:end-1,1:end-2).*dY)./(dY + dYW);
CzN = (Cz(2:end-1,2:end-1).*dZN + Cz(1:end-2,2:end-1).*dZ)./(dZ + dZN);
CzS = (Cz(2:end-1,2:end-1).*dZS + Cz(3:end-0,2:end-1).*dZ)./(dZ + dZS);

DE_m = (dY + dYE)./(dY./Dy(2:end-1,2:end-1) + dYE./Dy(2:end-1,3:end-0));
DW_m = (dY + dYW)./(dY./Dy(2:end-1,2:end-1) + dYW./Dy(2:end-1,1:end-2));
DS_m = (dZ + dZS)./(dZ./Dz(2:end-1,2:end-1) + dZS./Dz(3:end-0,2:end-1));
DN_m = (dZ + dZN)./(dZ./Dz(2:end-1,2:end-1) + dZN./Dz(1:end-2,2:end-1));

DE = [DE_m(:,1:end-1),Dy(2:end-1,end)];
DW = [Dy(2:end-1,1  ),DW_m(:,2:end-0)];
DN = [Dz(1  ,2:end-1);DN_m(2:end-0,:)];
DS = [DS_m(1:end-1,:);Dz(end,2:end-1)];

N = zeros(sz,sy);
S = zeros(sz,sy);
E = zeros(sz,sy);
W = zeros(sz,sy);
P = zeros(sz,sy);
R = zeros(sz,sy);

W = W + 2*dZ.*DW./(dY + dYW) - dZ.*CyW.*dY ./(dY + dYW);
P = P - 2*dZ.*DW./(dY + dYW) - dZ.*CyW.*dYW./(dY + dYW);

E = E + 2*dZ.*DE./(dY + dYE) + dZ.*CyE.*dY ./(dY + dYE);
P = P - 2*dZ.*DE./(dY + dYE) + dZ.*CyE.*dYE./(dY + dYE);

N = N + 2*dY.*DN./(dZ + dZN) - dY.*CzN.*dZ ./(dZ + dZN);
P = P - 2*dY.*DN./(dZ + dZN) - dY.*CzN.*dZN./(dZ + dZN);

S = S + 2*dY.*DS./(dZ + dZS) + dY.*CzS.*dZ ./(dZ + dZS);
P = P - 2*dY.*DS./(dZ + dZS) + dY.*CzS.*dZS./(dZ + dZS);

N(1  ,:) = 0;
W(:,1  ) = 0;
S(end,:) = 0;
E(:,end) = 0;

R(:,end) = R(:,end) - dZ(:,end).*(CyE(:,end).*dY(:,end)./(EB(1,:).'.*dY(:,end) + 2*EB(2,:).')).*EB(3,:).';
R(:,1  ) = R(:,1  ) + dZ(:,1  ).*(CyW(:,1  ).*dY(:,1  )./(WB(1,:).'.*dY(:,1  ) + 2*WB(2,:).')).*WB(3,:).';
R(end,:) = R(end,:) - dY(end,:).*(CzS(end,:).*dZ(end,:)./(SB(1,:)  .*dZ(end,:) + 2*SB(2,:)  )).*SB(3,:)  ;
R(1  ,:) = R(1  ,:) + dY(1  ,:).*(CzN(1  ,:).*dZ(1  ,:)./(NB(1,:)  .*dZ(1  ,:) + 2*NB(2,:)  )).*NB(3,:)  ;

R(:,end) = R(:,end) - dZ(:,end).*(2*DE(:,end)./(EB(1,:).'.*dY(:,end) + 2*EB(2,:).')).*EB(3,:).';
R(:,1  ) = R(:,1  ) - dZ(:,1  ).*(2*DW(:,1  )./(WB(1,:).'.*dY(:,1  ) + 2*WB(2,:).')).*WB(3,:).';
R(end,:) = R(end,:) - dY(end,:).*(2*DS(end,:)./(SB(1,:)  .*dZ(end,:) + 2*SB(2,:)  )).*SB(3,:)  ;
R(1  ,:) = R(1  ,:) - dY(1  ,:).*(2*DN(1  ,:)./(NB(1,:)  .*dZ(1  ,:) + 2*NB(2,:)  )).*NB(3,:)  ;

P(:,end) = P(:,end) + dZ(:,end).*(2*EB(2,:).'.*CyE(:,end))./(EB(1,:).'.*dY(:,end) + 2*EB(2,:).');
P(:,1  ) = P(:,1  ) - dZ(:,1  ).*(2*WB(2,:).'.*CyW(:,1  ))./(WB(1,:).'.*dY(:,1  ) + 2*WB(2,:).');
P(end,:) = P(end,:) + dY(end,:).*(2*SB(2,:)  .*CzS(end,:))./(SB(1,:)  .*dZ(end,:) + 2*SB(2,:)  );
P(1  ,:) = P(1  ,:) - dY(1  ,:).*(2*NB(2,:)  .*CzN(1  ,:))./(NB(1,:)  .*dZ(1  ,:) + 2*NB(2,:)  );

P(:,end) = P(:,end) + 2*dZ(:,end).*DE(:,end).*(2*EB(2,:).'./dY(:,end))./(EB(1,:).'.*dY(:,end) + 2*EB(2,:).');
P(:,1  ) = P(:,1  ) + 2*dZ(:,1  ).*DW(:,1  ).*(2*WB(2,:).'./dY(:,1  ))./(WB(1,:).'.*dY(:,1  ) + 2*WB(2,:).');
P(end,:) = P(end,:) + 2*dY(end,:).*DS(end,:).*(2*SB(2,:)  ./dZ(end,:))./(SB(1,:)  .*dZ(end,:) + 2*SB(2,:)  );
P(1  ,:) = P(1  ,:) + 2*dY(1  ,:).*DN(1  ,:).*(2*NB(2,:)  ./dZ(1  ,:))./(NB(1,:)  .*dZ(1  ,:) + 2*NB(2,:)  );

N = N./(dY.*dZ);
E = E./(dY.*dZ);
W = W./(dY.*dZ);
S = S./(dY.*dZ);
P = P./(dY.*dZ) + Prp(2:end-1,2:end-1);
R = R./(dY.*dZ) + RHS(2:end-1,2:end-1);

[s_z,s_y] = size(E);

e = reshape(E.',[1,s_y*s_z]);
w = reshape(W.',[1,s_y*s_z]);
n = reshape(N.',[1,s_y*s_z]);
s = reshape(S.',[1,s_y*s_z]);
p = reshape(P.',[1,s_y*s_z]);
r = reshape(R.',[1,s_y*s_z]);

e = e(1:end-1);
w = w(2:end-0);
n = n(s_y+1:end);
s = s(1:s_y*(s_z-1));

rows_p = 1:1:s_y*s_z;
rows_e = 1:1:s_y*s_z-1;
rows_w = 2:1:s_y*s_z;
rows_n = s_y+1:1:s_y*s_z;
rows_s = 1:1:s_y*s_z-s_y;

cols_p = 1:1:s_y*s_z;
cols_e = 2:1:s_y*s_z;
cols_w = 1:1:s_y*s_z-1;
cols_n = 1:1:s_y*s_z-s_y;
cols_s = s_y+1:1:s_y*s_z;

rows_vals = [rows_p, rows_e, rows_w, rows_n, rows_s];
cols_vals = [cols_p, cols_e, cols_w, cols_n, cols_s];
values = [p, e, w, n, s];

Coeffs = sparse(rows_vals, cols_vals, values);
Psi = (Coeffs\r.').';
U = reshape(Psi,[s_y,s_z]).';

N_val = (NB(3,:)  + 2*U(1  ,:).*NB(2,:) ./dY(1  ,:))./(NB(1,:)  + 2*NB(2,:) ./dY(1  ,:));
W_val = (WB(3,:)' + 2*U(:,1  ).*WB(2,:)'./dY(:,1  ))./(WB(1,:)' + 2*WB(2,:)'./dY(:,1  ));

S_val = (SB(3,:)  + 2*U(end,:).*SB(2,:) ./dZ(end,:))./(SB(1,:)  + 2*SB(2,:) ./dZ(end,:));
E_val = (EB(3,:)' + 2*U(:,end).*EB(2,:)'./dZ(:,end))./(EB(1,:)' + 2*EB(2,:)'./dZ(:,end));

U_ext = zeros(size(U)+2);

U_ext(2:end-1,2:end-1) = U;

U_ext(1  ,2:end-1) = N_val;
U_ext(2:end-1,1  ) = W_val;
U_ext(end,2:end-1) = S_val;
U_ext(2:end-1,end) = E_val;

U_ext(1,1) = N_val(1)/2+W_val(1)/2;

U_ext(end,1) = S_val(1)/2+W_val(end)/2;
U_ext(1,end) = N_val(end)/2+E_val(1)/2;

U_ext(end,end) = S_val(end)/2+E_val(end)/2;

end