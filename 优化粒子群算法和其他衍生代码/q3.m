clc
clear
T=200;
bbest=1;
X=[];
n=100 ; %粒子数目
v=0.01*ones(100,1);
w=0.4;
c1=1;
c2=2;
c3=3;
vmax=0.1;
b=0.01*(normrnd(0,2,100, 100)+50);%初始化粒子位置
pbest=min(b,[],2);
gbest=min(pbest);
dT=50*ones(100,1);

for t=1:T
    
    
    for i=1:100
         [dt,u]=Findbest(b(i,:));
        if u>0
            pbest(i)=u;
            dT(i)=dt;
            
        end
    end
    tbest=min(dT);
    gbestn=min(pbest(find(dT==tbest)))
    
    %初始化结束
    v=w.*v+(c1*rand(100,1).*(-b+pbest)/6+c2*rand(100,1).*(gbest-b)/6)*tbest/20;
    %更新速度
    v=0.8*v;
    for i=1:100
        if v(i)>vmax
            v(i)=vmax;
        end
    end
    b=b+v;
    
    
    if gbestn<gbest
        gbest=gbestn;
    end
 end
 


function [dt,bbest]=Findbest(b)
bbest=[];
dt=10;
    for i =1:size(b)
        fake(i)=1;
        t(i)=calculateb_t(b(i));
        [X,Y,F,qx,qy,theta]=calculate_t(b(i),t(i)); 
        in0=drawall(t(i),qx,qy,theta);
       
        [in1,t1]=divid(t(i)-10,t(i),b(i))
        
                
                if dt>t(i)-t1 
                    bbest=b(i)
                    dt=t(i)-t1;
                
                end
               
                
            
           
        
        
    end
    
   
    
end
%二分法
function [inm,tm]=divid(t0,t1,b)
    a=b/(2*pi);
    while t1-t0>0.001
    tm=(t1+t0)/2;
    [X0,Y0,F0,qx0,qy0,theta0]=calculate_t(b,t0);
    [X1,Y1,F1,qx1,qy1,theta1]=calculate_t(b,t1);
    [Xm,Ym,Fm,qxm,qym,thetam]=calculate_t(b,tm);
    in0=drawall(t0,qx0,qy0,theta0);
    in1=drawall(t1,qx1,qy1,theta1);
    inm=drawall(tm,qxm,qym,thetam);
    
    
    if inm==1
        break
        
        
    elseif inm==0
        t0=tm;
        
    end
    
    
end

end
function [X,Y,F,qx,qy,theta]=calculate_t(b,t) 

        anset=[];
        a=b/(2*pi);
        syms ano r0
        anset0=[];
        t=t/2;
        ans=vpasolve(r0^4-17.6*(r0^3)+(8.8^2+a^2)*r0^2-17.6*r0*a^2+(8.8^2)*(a^2)-(a^2)*(t^2)==0,r0);
        S_vpa = double(ans);
        y=real(S_vpa);
        anset=[anset,2*y(3)-8.8];
        
        %%第一个 L=286
       
        x =anset;
        ans1=vpasolve(x^2+ano^2-2*x*ano*cos(x/a-ano/a)==2.86^2,ano,[x,x+0.275]);
        S_vpa1 = double(ans1);
        anset1=S_vpa1;
        anset=[anset;S_vpa1];
        %%其他 L=165
        for k=0:52
            x=anset1;
                
            ans1=vpasolve(x^2+ano^2-2*x*ano*cos(x/a-ano/a)==1.65^2,ano,[x,x+0.275]);
            S_vpa1 = double(ans1);
            
            anset1=S_vpa1;
            anset=[anset;S_vpa1];
            
        end
        X=anset.*cos(anset/a);
        Y=anset.*sin(anset/a);
        F=[];
        for i =1:53
            F=[F;X(i,:)];
            F=[F;Y(i,:)];
        end
        qx=(X(1:52,:)+X(2:53,:))/2;
        qy=(Y(1:52,:)+Y(2:53,:))/2;
        theta=-atan((Y(2:53,:)-Y(1:52,:))./(X(2:53,:)-X(1:52,:)));
        dist=[];
        for i=1:52
        d=sqrt((X(i,:)-X(i+1,:)).^2+(Y(i,:)-Y(i+1,:)).^2);
        dist=[dist;d];
end 
end
function k=istbest(t,b)
    k=1;
    for i =0.1:1:10
        
        [d,qx0,qy0,theta0]=calculate_t(b);
        if d==0
            k=0;
            break
        end
        
        in=drawall(t-i,qx0,qy0,theta0);
        
        if in==1
            k=0;
            break
        end
    end
end
function anset=calculateb_t(c)

    anset=[];
    syms ano t
    for b=c
        
        a=b/(2*pi);
        r0=(4.5+8.8)/2;
        anset0=[];
        ans=vpasolve(r0^4-17.6*(r0^3)+(8.8^2+a^2)*r0^2-17.6*r0*a^2+(8.8^2)*(a^2)-(a^2)*(t^2)==0,t,[0,inf]);
        S_vpa = double(ans);
        y=real(S_vpa);
        anset=[anset;2*y];
    end
    
end


    
    
    %% 画
    
    function [p,vertices_absolute]=draw_segment(x_joint, y_joint, theta, L, W) 
        
        
        % 绘制一段身体或蛇头
        half_length = L/ 2;
        half_width = W/ 2;
       
        % 偏移值（骨骼相对于前端的偏移）
       
       
        
        % 旋转矩形
        R = [cos(theta), -sin(theta); sin(theta), cos(theta)];
        
        % 矩形顶点相对于中心的坐标
        vertices_relative = [-half_length, -half_width;
                             half_length, -half_width;
                             half_length,  half_width;
                            -half_length,  half_width];
                        
        % 旋转后顶点的绝对坐标
        
        vertices_absolute = (  vertices_relative * R) + [x_joint, y_joint];
        p=vertices_absolute(2,:);
        
        % 绘制矩形
        fill(vertices_absolute(:, 1), vertices_absolute(:, 2), 'r');
        hold on
    end

    %% 画龙，同时判断有没有撞
    function flag=drawall(t,qx,qy,theta)
    title(t)
    % 头位置和方向
        
    [x_head, y_head, theta_head] =deal(qx(1),qy(1),theta(1));
        
        % 绘制头
        [~,p]=draw_segment(x_head, y_head, theta_head, 3.41,0.3);
        
        % 更新每个身体段的位置和方向
        
        for i = 1:51
            
            [x_body, y_body, theta_body] = deal(qx(i),qy(i),theta(i));
            plot(x_body,y_body,'b*')
            [~,po]=draw_segment(x_body, y_body, theta_body, 2.2, 0.3);
            hold on
            if i>2
                flag=iscoll(p,po);
                if flag==1
                    break
                end
            end
            
         end
        
        % 刷新图像
        pause(0.0001)
        
        cla;
    end
    %% 判断函数
    function flag=iscoll(point,rect)
        [in,on]=inpolygon(point(:,1),point(:,2),[rect(:,1);rect(1,1)],[rect(:,2);rect(1,2)]);
        plot([rect(:,1);rect(1,1)],[rect(:,2);rect(1,2)])
        hold on 
        plot(point(1),point(2),'r*')
     
        flag=sum(in)+sum(on)
    end
