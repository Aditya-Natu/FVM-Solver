%Code developed by Aditya Natu%

function [Y, Z, dY, dZ, Y_ext, Z_ext] = mk_grid(y,z)

[~,sy] = size(y);
[~,sz] = size(z);

Y = repmat((y(1:end-1)+y(2:end)) /2,[sz-1,1]);
Z = repmat((z(1:end-1)+z(2:end))'/2,[1,sy-1]);

dY = repmat(y(2:end) - y(1:end-1), [sz-1,1]);
dZ = repmat(z(2:end)'-z(1:end-1)', [1,sy-1]);

Y_ext = repmat([y(1),(y(1:end-1)+y(2:end)) /2, y(end)],[sz+1,1]);
Z_ext = repmat([z(1);(z(1:end-1)+z(2:end))'/2; z(end)],[1,sy+1]);

end
