function Soln = PentDiagSol(N,E,W,S,P,R)

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
Soln = reshape(Psi,[s_y,s_z]).';

end