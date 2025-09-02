clear
%% 绘图仿真
a=0.55/(2*pi);
load H.mat

X=H.*cos(H/a);
Y=H.*sin(H/a);
F=[];
for i =1:225
F=[F;X(i,:)];
F=[F;Y(i,:)];
end



plot([0:300],X(1,1:301),'color','#1072BD','Linewidth', 1.5)

hold on
plot([170:300],X(51,1:301),'color','#77AE43','Linewidth', 1.5)
hold on

plot([290:300],X(101,1:301),'Color','#EDB021','Linewidth', 1.5)
hold on


%plot([0:300],X(151,1:301),'r-','Linewidth', 1.5)

%hold on
%plot([0:300],X(201,1:301),'Color','#7F318D','Linewidth', 1.5)

%hold on

