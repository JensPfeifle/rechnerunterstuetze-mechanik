function f = func(t,a)

if (a==1)
    f = cos(100*t);
elseif (a==2)
    f = t.*t;
elseif (a==3)
    f = exp(-5*t);
else
    f = 0;
end
