function err = ErrFind(Psi,Psi_old)

DiffMat = abs(Psi - Psi_old);

[hvals, rownum] = max(DiffMat);
[MaxDiff,colnum] = max(hvals);

err = abs(MaxDiff/Psi(rownum(colnum),colnum));

end