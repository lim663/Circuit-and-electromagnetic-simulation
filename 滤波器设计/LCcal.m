Z0=50;
fc=10^9;
gk=[1.1812,1.4228,2.0967,1.5734,2.0967,1.4228,1.1812];
res=[];
for i=1:7
    if rem(i,2)==0
        res(1,i)=Z0*gk(1,i)/(2*pi*fc);
    else
        res(1,i)=gk(1,i)/(2*pi*fc*Z0);
    end
end
