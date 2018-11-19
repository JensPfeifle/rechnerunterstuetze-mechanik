function df = dfunc(t,a)

if (a==1)
    df = -100*sin(100*t);
elseif (a==2)
    df = 2*t;
elseif (a==3)
    df = -5*exp(-5*t);
else
    df = 0;
end
