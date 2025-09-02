clc
clear
a=1.7/(2*pi);
%% 龙头
syms r0
%ans=solve(r0^4-17.6*(r0^3)+(8.8^2+a^2)*r0^2-17.6*r0*a^2+(8.8^2)*(a^2)-(a^2)*(t^2)==0,r0)
%S_vpa = double(ans)
anset=[];
anset0=[];
for t=0:0.5:444
    ans=solve(r0^4-17.6*(r0^3)+(8.8^2+a^2)*r0^2-17.6*r0*a^2+(8.8^2)*(a^2)-(a^2)*(t^2)==0,r0);
    S_vpa = double(ans);
    y=real(S_vpa);
    anset=[anset,y(3)];
    if mod(t*10,2)==0
        plot(t,2*anset(t+1)-8.8,'.','Color','r');
        anset0=[anset0,2*anset(t+1)-8.8];
        hold on 
    end
end
k=anset';
load H.mat
%% 迭代求其他板凳
Ben=[];%初始化
%%第一个 L=286
syms ano
anset1=[];

for x =anset0(300:445)
    
    ans1=vpasolve(x^2+ano^2-2*x*ano*cos(x/a-ano/a)==2.86^2,ano,[x,x+0.275]);
    S_vpa1 = double(ans1);
    anset1=[anset1,S_vpa1];

end



%%其他 L=165
for k=0:222
    anset2=[];
    for x=anset1
        
        ans1=vpasolve(x^2+ano^2-2*x*ano*cos(x/a-ano/a)==1.65^2,ano,[x,x+0.275]);
        S_vpa1 = double(ans1);
        anset2=[anset2,S_vpa1];
    end
    anset1=anset2;
    Ben=[Ben;anset2];
    plot(304,anset1)
end
load H.mat
%% 辛普森求微分算速度

DA=[];
R=sqrt(1+H.^2/a^2);
for i=3:431
    K=(-H(:,i+2)+8.*H(:,i+1)-8.*H(:,i-1)+H(:,i-2))/12;
    DA=[DA ,K];
end
K=(H(:,3)-H(:,1))/2;
DA=[K,DA];
K=(H(:,2)-H(:,1));
DA=[K,DA];
DAE=-DA.*R(:,1:431);
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
%% 绘图仿真


plot(X(:,1),Y(:,1),'g*-')

hold on
plot(X(:,60),Y(:,60),'r*-')

hold on
plot(X(:,120),Y(:,120),'y*-')

hold on
plot(X(:,180),Y(:,180),'p*-')

hold on
plot(X(:,400),Y(:,400),'b*-',LineWidth=70)
    
for t=1:400
    plot(X(:,t),Y(:,t),'r*-',LineWidth=7)
    pause(0.1)
    hold on
end

n=1:225
