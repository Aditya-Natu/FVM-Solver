%Developed and tested by Aditya Natu%

    function sol = ci_walls(a,b,n,r)

    if r == 1
        sol = linspace(a,b,n);
        
    else
        
        if mod(n,2)==1
            sol1 = geometric(a,(a+b)/2,(n+1)/2,r);   
            sol2 = geometric((a+b)/2,b,(n+1)/2,1/r);
            sol = [sol1(1:end-1),sol2];  
        else
            T1 = ((b-a)/2)/((r^(n/2 - 1)-1)/(r-1) + 0.5*r^(n/2 - 1));
            M = 0:1:n/2;
            d = T1*(r.^M - 1)/(r-1);
            c = T1*(r^(n/2) - 1)/(r-1);
            p = n/2-2:-1:0;
            e = T1*(r^(n/2 - 1) - r.^(p))/(r-1);
            sol = [a+d,a+c+e];
        end 
        
    end
    
    function sol = geometric(a,b,n,r)
        t1 = (b-a).*(1-r)./(1 - r.^(n-1));
        m = 0:1:n-1;
        sol = a + t1.*(1-r.^m)./(1-r);
    end
        
    end
