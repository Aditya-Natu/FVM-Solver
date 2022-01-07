function out = secant(f,x0,x1,tol,maxitn)

y0 = f(x0);
y1 = f(x1);

error = 1;
i=0;
while and(error>tol,i<maxitn)
    
    x = (x0*y1 - x1*y0)/(y1-y0);
    y = f(x);
    
    error = abs((x-x1)/x);

    x0 = x1;
    x1 = x;
    y0 = y1;
    y1 = y;
    
    i = i+1;
    
    fprintf('Iteration no. %f \n',i);
end

fprintf('Took %d iterations \n', i);

out = x;
