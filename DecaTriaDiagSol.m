function Soln = DecaTriaDiagSol(NN,NP,EE,EP,WW,WP,SS,SP,NE,NW,SE,SW,PP,R)

[s_z,s_y] = size(R);

ee = reshape(EE.',[1,s_y*s_z]);
ep = reshape(EP.',[1,s_y*s_z]);
ww = reshape(WW.',[1,s_y*s_z]);
wp = reshape(WP.',[1,s_y*s_z]);
nn = reshape(NN.',[1,s_y*s_z]);
np = reshape(NP.',[1,s_y*s_z]);
ss = reshape(SS.',[1,s_y*s_z]);
sp = reshape(SP.',[1,s_y*s_z]);
ne = reshape(NE.',[1,s_y*s_z]);
nw = reshape(NW.',[1,s_y*s_z]);
se = reshape(SE.',[1,s_y*s_z]);
sw = reshape(SW.',[1,s_y*s_z]);
pp = reshape(PP.',[1,s_y*s_z]);
r = reshape(R.',[1,s_y*s_z]);

ep = ep(1:end-1);
ee = ee(1:end-2);
wp = wp(2:end-0);
ww = ww(3:end-0);
np = np(1*s_y+1:end);
nn = nn(2*s_y+1:end);
sp = sp(1:s_y*(1*s_z-1));
ss = ss(1:s_y*(s_z-2));
ne = ne(1*s_y+1:end-1);
se = se(1:s_y*(1*s_z-1)-1);
nw = nw(1*s_y+2:end);
sw = sw(2:s_y*(1*s_z-1));

rows_pp = 1:1:s_y*s_z;
rows_ep = 1:1:s_y*s_z-1;
rows_ee = 1:1:s_y*s_z-2;
rows_wp = 2:1:s_y*s_z;
rows_ww = 3:1:s_y*s_z;
rows_np = 1*s_y+1:1:s_y*s_z;
rows_nn = 2*s_y+1:1:s_y*s_z;
rows_sp = 1:1:s_y*s_z-1*s_y;
rows_ss = 1:1:s_y*s_z-2*s_y;
rows_ne = 1*s_y+1:1:s_y*s_z-1;
rows_nw = 1*s_y+2:1:s_y*s_z;
rows_se = 1:1:s_y*s_z-1*s_y-1;
rows_sw = 2:1:s_y*s_z-1*s_y;


cols_pp = 1:1:s_y*s_z;
cols_ep = 2:1:s_y*s_z;
cols_ee = 3:1:s_y*s_z;
cols_wp = 1:1:s_y*s_z-1;
cols_ww = 1:1:s_y*s_z-2;
cols_np = 1:1:s_y*s_z-1*s_y;
cols_nn = 1:1:s_y*s_z-2*s_y;
cols_sp = 1*s_y+1:1:s_y*s_z;
cols_ss = 2*s_y+1:1:s_y*s_z;
cols_ne = 2:1:s_y*s_z-1*s_y;
cols_nw = 1:1:s_y*s_z-1*s_y-1;
cols_se = 1*s_y+2:1:s_y*s_z;
cols_sw = 1*s_y+1:1:s_y*s_z-1;

rows_vals = [rows_pp, rows_ep, rows_wp, rows_ee, rows_ww, rows_np, rows_ne, rows_nw, rows_sp, rows_se, rows_sw, rows_nn, rows_ss];
cols_vals = [cols_pp, cols_ep, cols_wp, cols_ee, cols_ww, cols_np, cols_ne, cols_nw, cols_sp, cols_se, cols_sw, cols_nn, cols_ss];
values = [pp, ep, wp, ee, ww, np, ne, nw, sp, se, sw, nn, ss];

Coeffs = sparse(rows_vals, cols_vals, values);
Psi = (Coeffs\r.').';
Soln = reshape(Psi,[s_y,s_z]).';

end