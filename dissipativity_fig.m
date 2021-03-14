clc
clear
%耗散output feedback控制 %non fragile controller
%%%模态转移阵的确定及其图形
finaltime=80;
p=[];
%Pi=[0.3 0.2 0.5;0.4 0.2 0.4;0.55 0.15 0.3];
p11=0.3;  p12=0.2;  p13=0.5;
p21=0.4;  p22=0.2;  p23=0.4;
p31=0.55; p32=0.15; p33=0.3; 

q=[];
Phi=[0.3 0.2 0.5;0.1 0.2 0.7;0.3 0.2 0.5];
q11=0.3; q12=0.2; q13=0.5;
q21=0.1; q22=0.2; q23=0.7;
q31=0.3; q32=0.2; q33=0.5;

flag  =1;                           %初始时刻（0时刻）的模态  
for i = 1:finaltime                 %1到500时刻的模态求取              
       if (flag == 1)
             a=rand;
             b=rand;
                if a<p11 
                    flag=1;           %system 1模态到1模态
                    if b<q11
                        flag1=1;
                    elseif q11<=b && b<=(q11+q12)
                        flag1=2;
                    else
                        flag1=3;
                    end
                                                
                 elseif  p11<=a && a<=(p11+p12)
                    flag=2;           %system 1模态到2模态
                    if b<q21
                        flag1=1;
                    elseif q21<=b && b<=(q21+q22)
                        flag1=2;
                    else
                        flag1=3;
                    end
                    
                else
                    flag=3;           %system 1模态到3模态 
                    if b<q31
                        flag1=1;
                    elseif q31<=b && b<=(q31+q32)
                        flag1=2;
                    else
                        flag1=3;
                    end
                            
                end
       elseif (flag==2)
             a=rand;
             b=rand;
                if a<p21
                    flag=1;           %system 2模态到1模态
                    if b<q11
                        flag1=1;
                    elseif q11<=b && b<=(q11+q12)
                        flag1=2;
                    else
                        flag1=3;
                    end
                elseif (p21<=a && a<=(p21+p22))
                    flag=2;           %system 2模态到2模态
                    if b<q21
                        flag1=1;
                    elseif q21<=b && b<=(q21+q22)
                        flag1=2;
                    else
                        flag1=3;
                    end

                else
                    flag=3;           %system2模态到3模态
                    if b<q31
                        flag1=1;
                    elseif q31<=b && b<=(q31+q32)
                        flag1=2;
                    else
                        flag1=3;
                    end
                end
       else
             a=rand;
             b=rand;
                if a<p31
                    flag=1;           %3模态到1模态
                    if b<q11
                        flag1=1;
                    elseif q11<=b && b<=(q11+q12)
                        flag1=2;
                    else
                        flag1=3;
                    end
                elseif (p31<=a && a<=(p31+p32))
                    flag=2;           %3模态到2模态
                    if b<q21
                        flag1=1;
                    elseif q21<=b && b<=(q21+q22)
                        flag1=2;
                    else
                        flag1=3;
                    end
                else
                    flag=3;           %3模态到3模态
                    if b<q31
                        flag1=1;
                    elseif q31<=b && b<=(q31+q32)
                        flag1=2;
                    else
                        flag1=3;
                    end
                end
             
       end
p(i)=flag;                              %1到500时刻的模态。一旦初始模态确定，后续的模态可以离线确定
q(i)=flag1;
end
p=[1 p];                                 %0到500时刻的模态。本质上，第500时刻的模态用不着
q=[1 q];
n1=0:finaltime;
figure
stairs(n1,p,'LineWidth',2,'Color','b');      %stairs函数绘制阶梯状图
%stem(n1,p,'fill','LineWidth',2,'Color','b')  %stem函数绘制火柴梗图
%stem(n1,q,'fill','LineWidth',1.5,'Color','k')%%黑色
axis([0 finaltime,0.5 3.5])
latexStr1 = ['$k$ '];
xlabel(latexStr1,'interpreter','latex','FontSize',14) 
latexStr2 = ['$\delta_k$ '];                 % LaTeX语法，因为有些东西不用latex写不出来
ylabel(latexStr2,'interpreter','latex','FontSize',14)

figure
stairs(n1,q,'LineWidth',2,'Color','b');      %stairs函数绘制阶梯状图
%stem(n1,p,'fill','LineWidth',2,'Color','b')  %stem函数绘制火柴梗图
%stem(n1,q,'fill','LineWidth',1.5,'Color','k')%%黑色
axis([0 finaltime,0.5 3.5])
latexStr1 = ['$k$ '];
xlabel(latexStr1,'interpreter','latex','FontSize',14) 
latexStr2 = ['$\eta_k$ '];                 % LaTeX语法，因为有些东西不用latex写不出来
ylabel(latexStr2,'interpreter','latex','FontSize',14)

%% system
%单连杆机械臂系统
M{1}=1; J{1}=1;    M{2}=5;J{2}=5;    M{3}=10;J{3}=10;    
L=0.5;g=9.81;beta1=0.01/pi;T=0.1;     %w取0.65时，tmin为10^-4级，不满足；w=1时，tmin为-10^-5级，满足
RR=2;        %区别QSR中的R

A{1,1}=[1 T;-T*g*L 1-T*RR/(J{1})];  A{2,1}=[1 T;-T*g*L 1-T*RR/(J{2})];  A{3,1}=[1 T;-T*g*L 1-T*RR/(J{3})];  
A{1,2}=[1 T;-beta1*T*g*L 1-T*RR/(J{1})];  A{2,2}=[1 T;-beta1*T*g*L 1-T*RR/(J{2})]; A{3,2}=[1 T;-beta1*T*g*L 1-T*RR/(J{3})]; 

B{1,1}=[0.1;T/(J{1})]; B{2,1}=[0;T/(J{2})];  B{3,1}=[0.5;T/(J{3})]; 
B{1,2}=[0;T/(J{1})];   B{2,2}=[0;T/(J{2})];  B{3,2}=[0.5;T/(J{3})]; 

D{1,1}=[0;T]; D{2,1}=[0;T]; D{3,1}=[0;T]; 
D{1,2}=[0;T]; D{2,2}=[0;T]; D{3,2}=[0;T]; 

C{1,1}=[1 0]; C{2,1}=[1 0]; C{3,1}=[1 0]; 
C{1,2}=[1 0]; C{2,2}=[1 0]; C{3,2}=[1 0]; 

E{1,1}=[1 0]; E{2,1}=[1 0]; E{3,1}=[1 0]; 
E{1,2}=[1 0]; E{2,2}=[1 0]; E{3,2}=[1 0]; 
%dissipativity
F{1,1}=[1]; F{2,1}=[1]; F{3,1}=[1]; 
F{1,2}=[1]; F{2,2}=[1]; F{3,2}=[1]; 


%l2-l-infinity
% F{1,1}=[0]; F{1,2}=[0]; F{1,3}=[0];  
% F{2,1}=[0]; F{2,2}=[0]; F{2,3}=[0]; 

%% 控制器增益
K{1,1}=-0.8644;  K{2,1}=-0.7739;  K{3,1}=-0.7278;  % dissipativity
K{1,2}=-1.7148;  K{2,2}=-1.6087;  K{3,2}=-1.5537;


% K{1}=-2.2727;  K{2}=-2.0759;  K{3}=-1.9698;  %l2-l-infinity
%% uncertainty delta_K=kesi1_vj*kesi2_vj(k)*kesi3_vj
kesi1{1,1}=0.02;  kesi1{2,1}=0.02;  kesi1{3,1}=0.02;
kesi1{1,2}=-0.02; kesi1{2,2}=-0.02; kesi1{3,2}=-0.02;

kesi3{1,1}=[0.02];  kesi3{2,1}=[0.02];  kesi3{3,1}=[0.02];
kesi3{1,2}=[-0.01]; kesi3{2,2}=[-0.01]; kesi3{3,2}=[-0.01];

%% 初值
x0=[0.5*pi;-0.1*pi];   %0时刻的初值
x(:,1)=x0;   %矩阵x中第一列元素是x0
x1(:,1)=x0;
y_hat0=[0.5*pi];    %0时刻估计的输出
y_hat(:,1)=y_hat0;
%%
ws=[];
for count=1:finaltime+1
    s=p(count);   %模态
    v=q(count);   %控制器模态  
    
    alpha1=[1 0]; prob = [0.8 0.2];   %prob 
    Alpha=randsrc(finaltime+1,1,[alpha1; prob]);   
    %randsrc(m,n,[alphabet;prob])参数m&n得到mxn的矩阵,alphabet是矩阵中元素可能取值的集合,prob是与alphabet同维度的概率向量，分别对应着alphabet中的每个元素可能出现的概率，其元素之和为1。
    beta=Alpha(count);  %beta(k)表示丢包的变量
    
    alpha=0.9;  %smoothing parameter
    
    x_1=[1 0]*x(:,count);     %x(:,count)表示第count列
    x_11=[1 0]*x1(:,count);  
            if x_1==0
                h_1=1;
            else
                h_1=(sin(x_1)-beta1*x_1)/(x_1*(1-beta1));
            end  
            %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
               if x_11==0
                 h_11=1;
             else
                 h_11=(sin(x_11)-beta1*x_11)/(x_11*(1-beta1));
             end
            h_2=1-h_1;
            h_22=1-h_11;
   
            %closed-loop system
   AA=h_1*A{s,1}+h_2*A{s,2}; 
   BB=h_1*B{s,1}+h_2*B{s,2};
   DD=h_1*D{s,1}+h_2*D{s,2};  
   
   CC=h_1*C{s,1}+h_2*C{s,2};
   EE=h_1*E{s,1}+h_2*E{s,2};
   FF=h_1*F{s,1}+h_2*F{s,2}; 
   
   %开环系统
   AA1=h_11*A{s,1}+h_22*A{s,2}; 
   DD1=h_11*D{s,1}+h_22*D{s,2};
   CC1=h_11*C{s,1}+h_22*C{s,2};
   EE1=h_11*E{s,1}+h_22*E{s,2};
   FF1=h_11*F{s,1}+h_22*F{s,2};
   
   %uncertainty
   kesi2_{1,1}=0.5*sin(count);   kesi2_{2,1}=0.5*sin(count);  kesi2_{3,1}=0.5*sin(count);
   kesi2_{1,2}=-0.5*cos(count);  kesi2_{2,2}=-0.5*cos(count); kesi2_{3,2}=-0.5*cos(count);
   
   Kh_{1}=h_1*(K{1,1}+kesi1{1,1}*kesi2_{1,1}*kesi3{1,1})+h_2*(K{1,2}+kesi1{1,2}*kesi2_{1,2}*kesi3{1,2});
   Kh_{2}=h_1*(K{2,1}+kesi1{2,1}*kesi2_{2,1}*kesi3{2,1})+h_2*(K{2,2}+kesi1{2,2}*kesi2_{2,2}*kesi3{2,2});
   Kh_{3}=h_1*(K{3,1}+kesi1{3,1}*kesi2_{3,1}*kesi3{3,1})+h_2*(K{3,2}+kesi1{3,2}*kesi2_{3,2}*kesi3{3,2});
  %没有控制器
    if count>=1 & count<=25
        w=0.5*exp(0.1*count)*sin(count);
        x1(:,count+1)=AA1*x1(:,count)+DD1*w;
        y1(:,count)=CC1*x1(:,count);
        z1(:,count)=EE1*x1(:,count)+FF1*w;
        
    else
        w=0;
        x1(:,count+1)=AA1*x1(:,count);
        y1(:,count)=CC1*x1(:,count);
        z1(:,count)=EE1*x1(:,count)+FF1*w;
    end  
 
%引入反馈控制器
    if count>=1 & count<=25
       w=0.5*exp(0.1*count)*sin(count);
       
     y_hat(:,count+1)=alpha*CC*x(:,count)+(1-alpha)*y_hat(:,count); %时序预测输出y_hat
     y(:,count)=CC*x(:,count);
     y_c(:,count)=beta*y(:,count)+(1-beta)*y_hat(:,count); %controller接收到的信号
     u(:,count)=Kh_{v}*y_c(:,count);
     x(:,count+1)=AA*x(:,count)+BB*u(:,count)+DD*w;
     z(:,count)=EE*x(:,count)+FF*w;
    else
       w=0;
       
     y_hat(:,count+1)=alpha*CC*x(:,count)+(1-alpha)*y_hat(:,count); %时序预测输出y_hat
     y(:,count)=EE*x(:,count);
     y_c(:,count)=beta*y(:,count)+(1-beta)*y_hat(:,count); %controller接收到的信号
     u(:,count)=Kh_{v}*y_c(:,count);
     x(:,count+1)=AA*x(:,count)+BB*u(:,count)+DD*w;
     z(:,count)=EE*x(:,count)+FF*w;
    
   end
        
end
figure    %开环系统的状态响应曲线 
%stairs(n2,q,'LineWidth',2,'Color','k');
plot([0:1:finaltime+1],x1,'LineWidth',2);
axis([0 finaltime,-4 4])
latexStr1 = ['$k$ '];
xlabel(latexStr1,'interpreter','latex','FontSize',14) 
latexStr2 = ['$State~trajectories$ '];
ylabel(latexStr2,'interpreter','latex','FontSize',14) 
latexStr1 = ['$x_{1}(k)$'];                 % LaTeX语法，因为有些东西不用latex写不出来
latexStr2 = ['$x_{2}(k)$'];
lgh=legend(latexStr1,latexStr2);  %加注释
set(lgh,'interpreter','latex','FontSize',16)    

figure           %闭环系统的状态响应曲线                            
plot([0:1:finaltime+1],x,'LineWidth',2);
axis([0 finaltime,-2 2])
latexStr1 = ['$k$ '];
xlabel(latexStr1,'interpreter','latex','FontSize',14) 
latexStr2 = ['$State~trajectories$ '];
ylabel(latexStr2,'interpreter','latex','FontSize',14)
%ylabel('State trajectories')
latexStr1 = ['$x_{1}(k)$'];                 % LaTeX语法，因为有些东西不用latex写不出来
latexStr2 = ['$x_{2}(k)$'];
lgh=legend(latexStr1,latexStr2);  %加注释
set(lgh,'interpreter','latex','FontSize',16)              %就是把字体改成latex格式figure

figure                                       
plot([0:1:finaltime],u,'LineWidth',2);
axis([0 finaltime,-4 2])
latexStr1 = ['$k$ '];
xlabel(latexStr1,'interpreter','latex','FontSize',14) 
latexStr2 = ['$u(k)$ '];                 % LaTeX语法，因为有些东西不用latex写不出来
ylabel(latexStr2,'interpreter','latex','FontSize',14)

% figure           %yuce曲线
% plot([0:1:finaltime],y,'LineWidth',2);
% % axis([0 finaltime,-1 1])
% latexStr1 = ['$k$ '];
% xlabel(latexStr1,'interpreter','latex','FontSize',14) 
% latexStr2 = ['$y(k)$ '];
% ylabel(latexStr2,'interpreter','latex','FontSize',14)
% % hold on           %prediction
% figure
% plot([0:1:finaltime+1],y_hat,'LineWidth',2);
% % axis([0 finaltime-1,-1 1])
% latexStr1 = ['$k$ '];                 % LaTeX语法，因为有些东西不用latex写不出来
% xlabel(latexStr1,'interpreter','latex','FontSize',14) 
% latexStr2 = ['$\hat{y}(k)$'];
% ylabel(latexStr2,'interpreter','latex','FontSize',14)


figure
plot([0:1:finaltime],y,'LineWidth',2);
hold on
plot([0:1:finaltime+1],y_hat,'LineWidth',2);
latexStr1 = ['$k$ '];
xlabel(latexStr1,'interpreter','latex','FontSize',14) 
latexStr2 = ['$outputs$ '];
ylabel(latexStr2,'interpreter','latex','FontSize',14) 
latexStr1 = ['$y(k)$'];                 % LaTeX语法，因为有些东西不用latex写不出来
latexStr2 = ['$\hat{y}(k)$'];
lgh=legend(latexStr1,latexStr2);  %加注释
set(lgh,'interpreter','latex','FontSize',16)    



