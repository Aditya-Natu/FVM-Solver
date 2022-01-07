function [Y, Z, dY, dZ, Y_ext, Z_ext, dY_ext, dZ_ext] = GridMaker(y,z)

[~,sy] = size(y);
[~,sz] = size(z);

Y = repmat((y(1:end-1)+y(2:end)) /2,[sz-1,1]);
Z = repmat((z(1:end-1)+z(2:end))'/2,[1,sy-1]);

dY = repmat(y(2:end) - y(1:end-1), [sz-1,1]);
dZ = repmat(z(2:end)'-z(1:end-1)', [1,sy-1]);

Y_ext = repmat([y(1),(y(1:end-1)+y(2:end)) /2, y(end)],[sz+1,1]);
Z_ext = repmat([z(1);(z(1:end-1)+z(2:end))'/2; z(end)],[1,sy+1]);

dY_int = [dY(1,:);dY;dY(end,:)];
dY_ext = [dY_int(:,1)/2,dY_int,dY_int(:,end)/2];

dZ_int = [dZ(:,1),dZ,dZ(:,end)];
dZ_ext = [dZ_int(1,:)/2;dZ_int;dZ_int(end,:)/2];

end