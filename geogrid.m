function sol = geogrid(a,b,n,r)

    sol1 = geometric(a,(a+b)/2,(n+1)/2,r);
    sol2 = geometric((a+b)/2,b,(n+1)/2,1/r);
    sol = [sol1(1:end-1),sol2]; 

end