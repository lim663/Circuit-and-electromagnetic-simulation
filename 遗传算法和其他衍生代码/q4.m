clc
clear
a=1.7/(2*pi);
%% 龙头
Ben=[];
%ans=solve(r0^4-17.6*(r0^3)+(4.5^2+a^2)*r0^2-17.6*r0*a^2+(4.5^2)*(a^2)-(a^2)*(t^2)==0,r0)
%S_vpa = double(ans)
anset=[];
anset0=[];
for t=-100:50:100
    
    if mod(t*10,2)==0
        plot(t,2*f(t/2)-4.5,'.','Color','r');
        anset=[anset,2*f(t/2)-4.5];

        hold on 
    end
end
Ben=[Ben;anset];
%初始化
%%第一个 L=286
syms ano
anset1=[];

for x =anset
    
    ans1=vpasolve(x^2+ano^2-2*x*ano*cos(x/a-ano/a)==2.86^2,ano);
    S_vpa1 = double(ans1(1));
    anset1=[anset1,S_vpa1];

end
Ben=[Ben;anset1];
k=anset';
for k=0:222
    anset2=[];
    for x=anset1
        
        ans1=vpasolve(x^2+ano^2-2*x*ano*cos(x/a-ano/a)==1.65^2,ano);
        S_vpa1 = double(ans1(1));
        anset2=[anset2,S_vpa1];
    end
    anset1=anset2;
    Ben=[Ben;anset2];
    
end
H=Ben;

X=H.*cos(H/a);
Y=H.*sin(H/a);
F=[];
for i =1:225
    F=[F;X(i,:)];
    F=[F;Y(i,:)];
end
%% 结果验证
dist=[];
for i=1:224
    d=sqrt((X(i,:)-X(i+1,:)).^2+(Y(i,:)-Y(i+1,:)).^2);
    dist=[dist;d];
end 
function S_vpa=f(t)
    syms r0
    a=1.7/(2*pi);
    ans=solve(r0^4-17.6*(r0^3)+(4.5^2+a^2)*r0^2-17.6*r0*a^2+(4.5^2)*(a^2)-(a^2)*(t^2)==0,r0);
    ans1 = real(ans);
    S_vpa=double(ans1(4));
    
end