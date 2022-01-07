function [U_ext] = BoundCond(NB,EB,WB,SB,dY,dZ,U)

N = (NB(3,:) + 2*U(1,:).*NB(2,:)./dY(1,:))./(NB(1,:) + 2*NB(2,:)./dY(1,:));
W = (WB(3,:)' + 2*U(:,1).*WB(2,:)'./dY(:,1))./(WB(1,:)' + 2*WB(2,:)'./dY(:,1));

S = (SB(3,:) + 2*U(end,:).*SB(2,:)./dZ(end,:))./(SB(1,:) + 2*SB(2,:)./dZ(end,:));
E = (EB(3,:)' + 2*U(:,end).*EB(2,:)'./dZ(:,end))./(EB(1,:)' + 2*EB(2,:)'./dZ(:,end));

U_ext = zeros(size(U)+2);

U_ext(2:end-1,2:end-1) = U;

U_ext(1,2:end-1) = N;
U_ext(2:end-1,1) = W;
U_ext(end,2:end-1) = S;
U_ext(2:end-1,end) = E;

U_ext(1,1) = N(1)/2+W(1)/2;

U_ext(end,1) = S(1)/2+W(end)/2;
U_ext(1,end) = N(end)/2+E(1)/2;

U_ext(end,end) = S(end)/2+E(end)/2;

end