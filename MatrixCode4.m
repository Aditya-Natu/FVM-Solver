function [NN,NP,EE,EP,WW,WP,SS,SP,NE,NW,SE,SW,PP,R] = MatrixCode4(dY, dZ)
 
 [sz,sy] = size(dY);
 
 dYPP = dY(2:end-1,2:end-1);
 dZPP = dZ(2:end-1,2:end-1);
 
 dYEP = dY(2:end-1,3:end  );
 dYEE = [dY(2:end-1,4:end),  zeros(sz-2,1)];
 dYWP = dY(2:end-1,1:end-2);
 dYWW = [zeros(sz-2,1),dY(2:end-1,1:end-3)];
 
 dZNP = dZ(1:end-2,2:end-1);
 dZNN = [zeros(1,sy-2);dZ(1:end-3,2:end-1)];
 dZSP = dZ(3:end  ,2:end-1);
 dZSS = [dZ(4:end,2:end-1);  zeros(1,sy-2)];
 
 ep = 2./(dYPP.*(dYPP + dYEP));
 wp = 2./(dYPP.*(dYPP + dYWP));
 np = 2./(dZPP.*(dZPP + dZNP));
 sp = 2./(dZPP.*(dZPP + dZSP));
 pp = - ep - wp - np - sp;
 
 ee = 2./(dYEP.*(dYEP + dYEE));
 we = 2./(dYEP.*(dYEP + dYPP));
 ne = 2./(dZPP.*(dZPP + dZNP));
 se = 2./(dZPP.*(dZPP + dZSP));
 pe = - ee - we - ne - se;
 
 ew = 2./(dYWP.*(dYWP + dYPP));
 ww = 2./(dYWP.*(dYWP + dYWW));
 nw = 2./(dZPP.*(dZPP + dZNP));
 sw = 2./(dZPP.*(dZPP + dZSP));
 pw = - ew - ww - nw - sw;
 
 en = 2./(dYPP.*(dYPP + dYEP));
 wn = 2./(dYPP.*(dYPP + dYWP));
 nn = 2./(dZNP.*(dZNP + dZNN));
 sn = 2./(dZNP.*(dZNP + dZPP));
 pn = - en - wn - nn - sn; 
 
 es = 2./(dYPP.*(dYPP + dYEP));
 ws = 2./(dYPP.*(dYPP + dYWP));
 ns = 2./(dZSP.*(dZSP + dZPP));
 ss = 2./(dZSP.*(dZSP + dZSS));
 ps = - es - ws - ns - ss;
 
 NP = np.*pp + pn.*np;
 SP = sp.*pp + ps.*sp;
 EP = ep.*pp + pe.*ep;
 WP = wp.*pp + pw.*wp;
 PP = pp.*pp + we.*ep + ew.*wp + ns.*sp + sn.*np;
 NN = nn.*np;
 SS = ss.*sp;
 EE = ee.*ep;
 WW = ww.*wp;
 NE = ne.*ep + en.*np;
 SE = se.*ep + es.*sp;
 NW = nw.*wp + wn.*np;
 SW = sw.*wp + ws.*sp;
 
 R = zeros(sz,sy);
 
 NP(1:2,:) = 0;
 SP(end-1:end,:) = 0;
 
 WP(:,1:2) = 0;
 EP(:,end-1:end) = 0;
 
 NN(1:2,:) = 0;
 SS(end-1:end,:) = 0;
 
 WW(:,1:2) = 0;
 EE(:,end-1:end) = 0;
 
 NE(1:2,:) = 0;
 NE(:,end-1:end) = 0;
 
 NW(1:2,:) = 0;
 NW(:,1:2) = 0;
 
 SE(end-1:end,:) = 0;
 SE(:,end-1:end) = 0;
 
 SW(end-1:end,:) = 0;
 SW(:,1:2) = 0;
 
end