%Lyapunov function V=x*P_sh*x
%dissipativity 
clear
clc
%fuzzy Markov jump system的dissipativity-based nonfragile输出反馈控制 
tic
format long
%% system
%单连杆机械臂系统
M{1}=1; J{1}=1;    M{2}=5;J{2}=5;    M{3}=10;J{3}=10;   
L=0.5;g=9.81;beta1=0.01/pi;T=0.1;     
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

F{1,1}=[1]; F{2,1}=[1]; F{3,1}=[1]; 
F{1,2}=[1]; F{2,2}=[1]; F{3,2}=[1]; 

%%%%%%%%%%%%%%%%%uncertainty%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
kesi1{1}=0.02;  kesi1{2}=0.02;  kesi1{3}=0.02;
 

 kesi3{1}=-0.02;  kesi3{2}=-0.02;  kesi3{3}=-0.02;
% kesi1{1}=0;  kesi1{2}=0;  kesi1{3}=0;
% kesi3{1}=0;  kesi3{2}=0;  kesi3{3}=0;

%% 转移矩阵 条件概率矩阵
rule=2;
S1=3;
S2=3;
Pi=[0.3 0.2 0.5;0.4 0.2 0.4;0.55 0.15 0.3];
% Phi=[0.3 0.2 0.5;0.1 0.2 0.7;0.3 0.2 0.5];
Phi=[0.5 0.2 0.3;0.1 0.3 0.6;0.4 0.1 0.5];

%% 耗散的参数
 Q=-1;    S=1;    R=1;  

Q1_=-chol(-Q,'lower');    % -Q=Q_'Q_
%% data dropouts
beta=0.7;
beta_=sqrt(beta*(1-beta));
%% 标量
aib1=0.1;
aib2=0.1;    %辅助求解LMIs
alpha=0.8;  %平滑系数
%% 变量
setlmis([]); 
for i=1:rule
    for s=1:S1 
        for n=1:S2
            P1(s,i)=lmivar(1,[2,1]);       %2x2对称矩阵
            P2(s,i)=lmivar(2,[2,1]);
            P3(s,i)=lmivar(1,[1,0]);
            
            W1(s,n,i)=lmivar(1,[2,1]);  
            W2(s,n,i)=lmivar(2,[2,1]);
            W3(s,n,i)=lmivar(1,[1,0]);
        end
    end 
end


for n=1:S2
    Y(n)=lmivar(1,[1,0]);         %1x1的矩阵
    X(n)=lmivar(1,[1,0]);
    
end

gamma=lmivar(1,[1,0]);           %标量   %耗散指标
%%  保证P>0  Q>0
for i=1:rule
     for s=1:S1
     FF(s,i)=newlmi;         %用函数newlmi确定线性不等式的名称为FF(m,i)
     lmiterm([FF(s,i) 1 1 P1(s,i)],-1,1);   
     lmiterm([FF(s,i) 1 2 P2(s,i)],-1,1); 
     lmiterm([FF(s,i) 2 2 P3(s,i)],-1,1); 
     end    
end

for i=1:rule
    for s=1:S1
        for n=1:S2
        FE(s,n,i)=newlmi;         
        lmiterm([FE(s,n,i) 1 1 W1(s,n,i)],-1,1);
        lmiterm([FE(s,n,i) 1 2 W2(s,n,i)],-1,1);
        lmiterm([FE(s,n,i) 2 2 W3(s,n,i)],-1,1);
        
        end       
    end
end

FEA=newlmi;
lmiterm([FEA 1 1 gamma],-1,1);
%% 定理2中的（28）
for i=1:rule
    for s=1:S1   
        LMI{s,i}=newlmi;
        lmiterm([LMI{s,i} 1 1 P1(s,i)],-1,1);
        lmiterm([LMI{s,i} 1 1 W1(s,1,i)],Phi(s,1),1);
        lmiterm([LMI{s,i} 1 1 W1(s,2,i)],Phi(s,2),1);                   
        lmiterm([LMI{s,i} 1 1 W1(s,3,i)],Phi(s,3),1); 
        
        lmiterm([LMI{s,i} 1 2 P2(s,i)],-1,1);
        lmiterm([LMI{s,i} 1 2 W2(s,1,i)],Phi(s,1),1);
        lmiterm([LMI{s,i} 1 2 W2(s,2,i)],Phi(s,2),1);                   
        lmiterm([LMI{s,i} 1 2 W2(s,3,i)],Phi(s,3),1); 
        
        lmiterm([LMI{s,i} 2 2 P3(s,i)],-1,1);
        lmiterm([LMI{s,i} 2 2 W3(s,1,i)],Phi(s,1),1);
        lmiterm([LMI{s,i} 2 2 W3(s,2,i)],Phi(s,2),1);                   
        lmiterm([LMI{s,i} 2 2 W3(s,3,i)],Phi(s,3),1); 
        
    end
end
% 定理2中的(29)式 ii
for f=1:rule
    for i=1:rule
        for s=1:S1 
            for v=1:S2
            LMI{f,i,s,v}=newlmi;
            lmiterm([LMI{f,i,s,v} 1 1 W1(s,v,i)],-1,1);
            lmiterm([LMI{f,i,s,v} 1 2 W2(s,v,i)],-1,1);
            lmiterm([LMI{f,i,s,v} 2 2 W3(s,v,i)],-1,1);
      
            lmiterm([LMI{f,i,s,v} 3 1 0],-S'*E{s,i});
            
            lmiterm([LMI{f,i,s,v} 3 3 0],-S'*F{s,i});
            lmiterm([LMI{f,i,s,v} 3 3 0],-F{s,i}'*S);
            lmiterm([LMI{f,i,s,v} 3 3 0],-R);
            lmiterm([LMI{f,i,s,v} 3 3 gamma],1,1);
            
            lmiterm([LMI{f,i,s,v} 4 1 0],Q1_*E{s,i});
            lmiterm([LMI{f,i,s,v} 4 3 0],Q1_*F{s,i});
            lmiterm([LMI{f,i,s,v} 4 4 0],-1);
            
            %omga41
            lmiterm([LMI{f,i,s,v} 5 1 P1(1,f)],1,sqrt(Pi(s,1))*A{s,i});
            lmiterm([LMI{f,i,s,v} 5 1 P2(1,f)],1,sqrt(Pi(s,1))*alpha*C{s,i});
            lmiterm([LMI{f,i,s,v} 5 2 P2(1,f)],1,sqrt(Pi(s,1))*(1-alpha));
            lmiterm([LMI{f,i,s,v} 6 1 -P2(1,f)],1,sqrt(Pi(s,1))*A{s,i});
            lmiterm([LMI{f,i,s,v} 6 1 P3(1,f)],1,sqrt(Pi(s,1))*alpha*C{s,i});
            lmiterm([LMI{f,i,s,v} 6 2 P3(1,f)],1,sqrt(Pi(s,1))*(1-alpha));           
            lmiterm([LMI{f,i,s,v} 5 1 Y(v)],sqrt(Pi(s,1))*beta*B{s,i},C{s,i});
            lmiterm([LMI{f,i,s,v} 5 2 Y(v)],sqrt(Pi(s,1))*(1-beta)*B{s,i},1);
            
            lmiterm([LMI{f,i,s,v} 7 1 P1(2,f)],1,sqrt(Pi(s,2))*A{s,i});
            lmiterm([LMI{f,i,s,v} 7 1 P2(2,f)],1,sqrt(Pi(s,2))*alpha*C{s,i});
            lmiterm([LMI{f,i,s,v} 7 2 P2(2,f)],1,sqrt(Pi(s,2))*(1-alpha));
            lmiterm([LMI{f,i,s,v} 8 1 -P2(2,f)],1,sqrt(Pi(s,2))*A{s,i});
            lmiterm([LMI{f,i,s,v} 8 1 P3(2,f)],1,sqrt(Pi(s,2))*alpha*C{s,i});
            lmiterm([LMI{f,i,s,v} 8 2 P3(2,f)],1,sqrt(Pi(s,2))*(1-alpha));
            lmiterm([LMI{f,i,s,v} 7 1 Y(v)],sqrt(Pi(s,2))*beta*B{s,i},C{s,i});
            lmiterm([LMI{f,i,s,v} 7 2 Y(v)],sqrt(Pi(s,2))*(1-beta)*B{s,i},1);
            
            lmiterm([LMI{f,i,s,v} 9 1 P1(3,f)],1,sqrt(Pi(s,3))*A{s,i});
            lmiterm([LMI{f,i,s,v} 9 1 P2(3,f)],1,sqrt(Pi(s,3))*alpha*C{s,i});
            lmiterm([LMI{f,i,s,v} 9 2 P2(3,f)],1,sqrt(Pi(s,3))*(1-alpha));
            lmiterm([LMI{f,i,s,v} 10 1 -P2(3,f)],1,sqrt(Pi(s,3))*A{s,i});
            lmiterm([LMI{f,i,s,v} 10 1 P3(3,f)],1,sqrt(Pi(s,3))*alpha*C{s,i});
            lmiterm([LMI{f,i,s,v} 10 2 P3(3,f)],1,sqrt(Pi(s,3))*(1-alpha));
            lmiterm([LMI{f,i,s,v} 9 1 Y(v)],sqrt(Pi(s,3))*beta*B{s,i},C{s,i});
            lmiterm([LMI{f,i,s,v} 9 2 Y(v)],sqrt(Pi(s,3))*(1-beta)*B{s,i},1);
            
            %omga42
            lmiterm([LMI{f,i,s,v} 5 3 P1(1,f)],1,sqrt(Pi(s,1))*D{s,i});
            lmiterm([LMI{f,i,s,v} 6 3 -P2(1,f)],1,sqrt(Pi(s,1))*D{s,i});
            
            lmiterm([LMI{f,i,s,v} 7 3 P1(2,f)],1,sqrt(Pi(s,2))*D{s,i});
            lmiterm([LMI{f,i,s,v} 8 3 -P2(2,f)],1,sqrt(Pi(s,2))*D{s,i});
            
            lmiterm([LMI{f,i,s,v} 9 3 P1(3,f)],1,sqrt(Pi(s,3))*D{s,i});
            lmiterm([LMI{f,i,s,v} 10 3 -P2(3,f)],1,sqrt(Pi(s,3))*D{s,i});
            
            %-omga44
            lmiterm([LMI{f,i,s,v} 5 5 P1(1,f)],-1,1);
            lmiterm([LMI{f,i,s,v} 5 6 P2(1,f)],-1,1);
            lmiterm([LMI{f,i,s,v} 6 6 P3(1,f)],-1,1);
            
            lmiterm([LMI{f,i,s,v} 7 7 P1(2,f)],-1,1);
            lmiterm([LMI{f,i,s,v} 7 8 P2(2,f)],-1,1);
            lmiterm([LMI{f,i,s,v} 8 8 P3(2,f)],-1,1);
            
            lmiterm([LMI{f,i,s,v} 9 9 P1(3,f)],-1,1);
            lmiterm([LMI{f,i,s,v} 9 10 P2(3,f)],-1,1);
            lmiterm([LMI{f,i,s,v} 10 10 P3(3,f)],-1,1);
            
            %omga51
            lmiterm([LMI{f,i,s,v} 11 1 Y(v)],sqrt(Pi(s,1))*beta_*B{s,i},C{s,i});
            lmiterm([LMI{f,i,s,v} 11 2 Y(v)],-sqrt(Pi(s,1))*beta_*B{s,i},1);
            
            lmiterm([LMI{f,i,s,v} 11 11 P1(1,f)],-1,1);
            lmiterm([LMI{f,i,s,v} 11 12 P2(1,f)],-1,1);
            lmiterm([LMI{f,i,s,v} 12 12 P3(1,f)],-1,1);
            
            lmiterm([LMI{f,i,s,v} 13 1 Y(v)],sqrt(Pi(s,2))*beta_*B{s,i},C{s,i});
            lmiterm([LMI{f,i,s,v} 13 2 Y(v)],-sqrt(Pi(s,2))*beta_*B{s,i},1);
            
            lmiterm([LMI{f,i,s,v} 13 13 P1(2,f)],-1,1);
            lmiterm([LMI{f,i,s,v} 13 14 P2(2,f)],-1,1);
            lmiterm([LMI{f,i,s,v} 14 14 P3(2,f)],-1,1);
            
            lmiterm([LMI{f,i,s,v} 15 1 Y(v)],sqrt(Pi(s,3))*beta_*B{s,i},C{s,i});
            lmiterm([LMI{f,i,s,v} 15 2 Y(v)],-sqrt(Pi(s,3))*beta_*B{s,i},1);
            
            lmiterm([LMI{f,i,s,v} 15 15 P1(3,f)],-1,1);
            lmiterm([LMI{f,i,s,v} 15 16 P2(3,f)],-1,1);
            lmiterm([LMI{f,i,s,v} 16 16 P3(3,f)],-1,1);
            
            
            lmiterm([LMI{f,i,s,v} 17 1 0],kesi3{v}*C{s,i});
            lmiterm([LMI{f,i,s,v} 17 17 0],-aib1);
            lmiterm([LMI{f,i,s,v} 18 2 0],kesi3{v});
            lmiterm([LMI{f,i,s,v} 18 18 0],-aib1);
            
            lmiterm([LMI{f,i,s,v} 19 5 P1(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*beta*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 19 6 P2(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*beta*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 20 5 P1(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*(1-beta)*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 20 6 P2(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*(1-beta)*B{s,i}',1);
            
            lmiterm([LMI{f,i,s,v} 19 7 P1(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*beta*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 19 8 P2(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*beta*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 20 7 P1(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*(1-beta)*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 20 8 P2(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*(1-beta)*B{s,i}',1);
            
            lmiterm([LMI{f,i,s,v} 19 9 P1(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*beta*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 19 10 P2(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*beta*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 20 9 P1(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*(1-beta)*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 20 10 P2(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*(1-beta)*B{s,i}',1);
            
            lmiterm([LMI{f,i,s,v} 19 11 P1(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*beta_*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 19 12 P2(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*beta_*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 20 11 P1(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*(-beta_)*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 20 12 P2(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*(-beta_)*B{s,i}',1);
            
            lmiterm([LMI{f,i,s,v} 19 13 P1(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*beta_*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 19 14 P2(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*beta_*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 20 13 P1(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*(-beta_)*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 20 14 P2(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*(-beta_)*B{s,i}',1);
            
            lmiterm([LMI{f,i,s,v} 19 15 P1(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*beta_*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 19 16 P2(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*beta_*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 20 15 P1(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*(-beta_)*B{s,i}',1);
            lmiterm([LMI{f,i,s,v} 20 16 P2(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*(-beta_)*B{s,i}',1);
            
            lmiterm([LMI{f,i,s,v} 19 19 0],-aib1);
            lmiterm([LMI{f,i,s,v} 20 20 0],-aib1);
            
            
            lmiterm([LMI{f,i,s,v} 21 1 Y(v)],aib2,C{s,i});
            lmiterm([LMI{f,i,s,v} 22 2 Y(v)],aib2,1);
            
            lmiterm([LMI{f,i,s,v} 21 5 P1(1,f)],sqrt(Pi(s,1))*beta*(B{s,i})',1);
            lmiterm([LMI{f,i,s,v} 21 5 X(v)],-sqrt(Pi(s,1))*beta,(B{s,i})');
            lmiterm([LMI{f,i,s,v} 21 6 P2(1,f)],sqrt(Pi(s,1))*beta*(B{s,i})',1);
            lmiterm([LMI{f,i,s,v} 22 5 P1(1,f)],sqrt(Pi(s,1))*(1-beta)*(B{s,i})',1);
            lmiterm([LMI{f,i,s,v} 22 5 X(v)],-sqrt(Pi(s,1))*(1-beta),(B{s,i})');
            lmiterm([LMI{f,i,s,v} 22 6 P2(1,f)],sqrt(Pi(s,1))*(1-beta)*(B{s,i})',1);
            
            lmiterm([LMI{f,i,s,v} 21 7 P1(2,f)],sqrt(Pi(s,2))*beta*(B{s,i})',1);
            lmiterm([LMI{f,i,s,v} 21 7 X(v)],-sqrt(Pi(s,2))*beta,(B{s,i})');
            lmiterm([LMI{f,i,s,v} 21 8 P2(2,f)],sqrt(Pi(s,2))*beta*(B{s,i})',1);
            lmiterm([LMI{f,i,s,v} 22 7 P1(2,f)],sqrt(Pi(s,2))*(1-beta)*(B{s,i})',1);
            lmiterm([LMI{f,i,s,v} 22 7 X(v)],-sqrt(Pi(s,2))*(1-beta),(B{s,i})');
            lmiterm([LMI{f,i,s,v} 22 8 P2(2,f)],sqrt(Pi(s,2))*(1-beta)*(B{s,i})',1);
            
            lmiterm([LMI{f,i,s,v} 21 9 P1(3,f)],sqrt(Pi(s,3))*beta*(B{s,i})',1);
            lmiterm([LMI{f,i,s,v} 21 9 X(v)],-sqrt(Pi(s,3))*beta,(B{s,i})');
            lmiterm([LMI{f,i,s,v} 21 10 P2(3,f)],sqrt(Pi(s,3))*beta*(B{s,i})',1);
            lmiterm([LMI{f,i,s,v} 22 9 P1(3,f)],sqrt(Pi(s,3))*(1-beta)*(B{s,i})',1);
            lmiterm([LMI{f,i,s,v} 22 9 X(v)],-sqrt(Pi(s,3))*(1-beta),(B{s,i})');
            lmiterm([LMI{f,i,s,v} 22 10 P2(3,f)],sqrt(Pi(s,3))*(1-beta)*(B{s,i})',1);
            
            %omga66
            lmiterm([LMI{f,i,s,v} 21 21 X(v)],aib2,-1,'s');
            lmiterm([LMI{f,i,s,v} 22 22 X(v)],aib2,-1,'s');
            
            end
        end
    end
end
% 定理2中的(30)式 j<i
for f=1:rule
    for j=1:(rule-1)
        for i=(j+1):rule 
            for s=1:S1 
                for v=1:S2
            LMI{f,i,j,s,v}=newlmi;
    %%%%%%%%%%%%%ij
            lmiterm([LMI{f,i,j,s,v} 1 1 W1(s,v,i)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 1 2 W2(s,v,i)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 2 2 W3(s,v,i)],-1,1);
      
            lmiterm([LMI{f,i,j,s,v} 3 1 0],-S'*E{s,i});
            
            lmiterm([LMI{f,i,j,s,v} 3 3 0],-S'*F{s,i});
            lmiterm([LMI{f,i,j,s,v} 3 3 0],-F{s,i}'*S);
            lmiterm([LMI{f,i,j,s,v} 3 3 0],-R);
            lmiterm([LMI{f,i,j,s,v} 3 3 gamma],1,1);
            
            lmiterm([LMI{f,i,j,s,v} 4 1 0],Q1_*E{s,i});
            lmiterm([LMI{f,i,j,s,v} 4 3 0],Q1_*F{s,i});
            lmiterm([LMI{f,i,j,s,v} 4 4 0],-1);
            
            %omga41
            lmiterm([LMI{f,i,j,s,v} 5 1 P1(1,f)],1,sqrt(Pi(s,1))*A{s,i});
            lmiterm([LMI{f,i,j,s,v} 5 1 P2(1,f)],1,sqrt(Pi(s,1))*alpha*C{s,j});
            lmiterm([LMI{f,i,j,s,v} 5 2 P2(1,f)],1,sqrt(Pi(s,1))*(1-alpha));
            lmiterm([LMI{f,i,j,s,v} 6 1 -P2(1,f)],1,sqrt(Pi(s,1))*A{s,i});
            lmiterm([LMI{f,i,j,s,v} 6 1 P3(1,f)],1,sqrt(Pi(s,1))*alpha*C{s,j});
            lmiterm([LMI{f,i,j,s,v} 6 2 P3(1,f)],1,sqrt(Pi(s,1))*(1-alpha));           
            lmiterm([LMI{f,i,j,s,v} 5 1 Y(v)],sqrt(Pi(s,1))*beta*B{s,i},C{s,j});
            lmiterm([LMI{f,i,j,s,v} 5 2 Y(v)],sqrt(Pi(s,1))*(1-beta)*B{s,i},1);
            
            lmiterm([LMI{f,i,j,s,v} 7 1 P1(2,f)],1,sqrt(Pi(s,2))*A{s,i});
            lmiterm([LMI{f,i,j,s,v} 7 1 P2(2,f)],1,sqrt(Pi(s,2))*alpha*C{s,j});
            lmiterm([LMI{f,i,j,s,v} 7 2 P2(2,f)],1,sqrt(Pi(s,2))*(1-alpha));
            lmiterm([LMI{f,i,j,s,v} 8 1 -P2(2,f)],1,sqrt(Pi(s,2))*A{s,i});
            lmiterm([LMI{f,i,j,s,v} 8 1 P3(2,f)],1,sqrt(Pi(s,2))*alpha*C{s,j});
            lmiterm([LMI{f,i,j,s,v} 8 2 P3(2,f)],1,sqrt(Pi(s,2))*(1-alpha));
            lmiterm([LMI{f,i,j,s,v} 7 1 Y(v)],sqrt(Pi(s,2))*beta*B{s,i},C{s,j});
            lmiterm([LMI{f,i,j,s,v} 7 2 Y(v)],sqrt(Pi(s,2))*(1-beta)*B{s,i},1);
            
            lmiterm([LMI{f,i,j,s,v} 9 1 P1(3,f)],1,sqrt(Pi(s,3))*A{s,i});
            lmiterm([LMI{f,i,j,s,v} 9 1 P2(3,f)],1,sqrt(Pi(s,3))*alpha*C{s,j});
            lmiterm([LMI{f,i,j,s,v} 9 2 P2(3,f)],1,sqrt(Pi(s,3))*(1-alpha));
            lmiterm([LMI{f,i,j,s,v} 10 1 -P2(3,f)],1,sqrt(Pi(s,3))*A{s,i});
            lmiterm([LMI{f,i,j,s,v} 10 1 P3(3,f)],1,sqrt(Pi(s,3))*alpha*C{s,j});
            lmiterm([LMI{f,i,j,s,v} 10 2 P3(3,f)],1,sqrt(Pi(s,3))*(1-alpha));
            lmiterm([LMI{f,i,j,s,v} 9 1 Y(v)],sqrt(Pi(s,3))*beta*B{s,i},C{s,j});
            lmiterm([LMI{f,i,j,s,v} 9 2 Y(v)],sqrt(Pi(s,3))*(1-beta)*B{s,i},1);
            
            %omga42
            lmiterm([LMI{f,i,j,s,v} 5 3 P1(1,f)],1,sqrt(Pi(s,1))*D{s,i});
            lmiterm([LMI{f,i,j,s,v} 6 3 -P2(1,f)],1,sqrt(Pi(s,1))*D{s,i});
            
            lmiterm([LMI{f,i,j,s,v} 7 3 P1(2,f)],1,sqrt(Pi(s,2))*D{s,i});
            lmiterm([LMI{f,i,j,s,v} 8 3 -P2(2,f)],1,sqrt(Pi(s,2))*D{s,i});
            
            lmiterm([LMI{f,i,j,s,v} 9 3 P1(3,f)],1,sqrt(Pi(s,3))*D{s,i});
            lmiterm([LMI{f,i,j,s,v} 10 3 -P2(3,f)],1,sqrt(Pi(s,3))*D{s,i});
            
            %-omga44
            lmiterm([LMI{f,i,j,s,v} 5 5 P1(1,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 5 6 P2(1,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 6 6 P3(1,f)],-1,1);
            
            lmiterm([LMI{f,i,j,s,v} 7 7 P1(2,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 7 8 P2(2,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 8 8 P3(2,f)],-1,1);
            
            lmiterm([LMI{f,i,j,s,v} 9 9 P1(3,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 9 10 P2(3,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 10 10 P3(3,f)],-1,1);
            
            %omga51
            lmiterm([LMI{f,i,j,s,v} 11 1 Y(v)],sqrt(Pi(s,1))*beta_*B{s,i},C{s,j});
            lmiterm([LMI{f,i,j,s,v} 11 2 Y(v)],-sqrt(Pi(s,1))*beta_*B{s,i},1);
            
            
            lmiterm([LMI{f,i,j,s,v} 11 11 P1(1,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 11 12 P2(1,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 12 12 P3(1,f)],-1,1);
            
            lmiterm([LMI{f,i,j,s,v} 13 1 Y(v)],sqrt(Pi(s,2))*beta_*B{s,i},C{s,j});
            lmiterm([LMI{f,i,j,s,v} 13 2 Y(v)],-sqrt(Pi(s,2))*beta_*B{s,i},1);
            
            lmiterm([LMI{f,i,j,s,v} 13 13 P1(2,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 13 14 P2(2,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 14 14 P3(2,f)],-1,1);
            
            lmiterm([LMI{f,i,j,s,v} 15 1 Y(v)],sqrt(Pi(s,3))*beta_*B{s,i},C{s,j});
            lmiterm([LMI{f,i,j,s,v} 15 2 Y(v)],-sqrt(Pi(s,3))*beta_*B{s,i},1);
            
            lmiterm([LMI{f,i,j,s,v} 15 15 P1(3,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 15 16 P2(3,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 16 16 P3(3,f)],-1,1);
            
            lmiterm([LMI{f,i,j,s,v} 17 1 0],kesi3{v}*C{s,j});
            lmiterm([LMI{f,i,j,s,v} 17 17 0],-aib1);
            lmiterm([LMI{f,i,j,s,v} 18 2 0],kesi3{v});
            lmiterm([LMI{f,i,j,s,v} 18 18 0],-aib1);
            
            lmiterm([LMI{f,i,j,s,v} 19 5 P1(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*beta*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 19 6 P2(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*beta*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 20 5 P1(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*(1-beta)*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 20 6 P2(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*(1-beta)*B{s,i}',1);
            
            lmiterm([LMI{f,i,j,s,v} 19 7 P1(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*beta*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 19 8 P2(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*beta*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 20 7 P1(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*(1-beta)*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 20 8 P2(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*(1-beta)*B{s,i}',1);
            
            lmiterm([LMI{f,i,j,s,v} 19 9 P1(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*beta*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 19 10 P2(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*beta*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 20 9 P1(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*(1-beta)*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 20 10 P2(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*(1-beta)*B{s,i}',1);
            
            lmiterm([LMI{f,i,j,s,v} 19 11 P1(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*beta_*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 19 12 P2(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*beta_*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 20 11 P1(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*(-beta_)*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 20 12 P2(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*(-beta_)*B{s,i}',1);
            
            lmiterm([LMI{f,i,j,s,v} 19 13 P1(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*beta_*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 19 14 P2(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*beta_*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 20 13 P1(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*(-beta_)*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 20 14 P2(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*(-beta_)*B{s,i}',1);
            
            lmiterm([LMI{f,i,j,s,v} 19 15 P1(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*beta_*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 19 16 P2(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*beta_*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 20 15 P1(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*(-beta_)*B{s,i}',1);
            lmiterm([LMI{f,i,j,s,v} 20 16 P2(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*(-beta_)*B{s,i}',1);
            
            lmiterm([LMI{f,i,j,s,v} 19 19 0],-aib1);
            lmiterm([LMI{f,i,j,s,v} 20 20 0],-aib1);
            
            
            lmiterm([LMI{f,i,j,s,v} 21 1 Y(v)],aib2,C{s,j});
            lmiterm([LMI{f,i,j,s,v} 22 2 Y(v)],aib2,1);
            
            lmiterm([LMI{f,i,j,s,v} 21 5 P1(1,f)],sqrt(Pi(s,1))*beta*(B{s,i})',1);
            lmiterm([LMI{f,i,j,s,v} 21 5 X(v)],-sqrt(Pi(s,1))*beta,(B{s,i})');
            lmiterm([LMI{f,i,j,s,v} 21 6 P2(1,f)],sqrt(Pi(s,1))*beta*(B{s,i})',1);
            lmiterm([LMI{f,i,j,s,v} 22 5 P1(1,f)],sqrt(Pi(s,1))*(1-beta)*(B{s,i})',1);
            lmiterm([LMI{f,i,j,s,v} 22 5 X(v)],-sqrt(Pi(s,1))*(1-beta),(B{s,i})');
            lmiterm([LMI{f,i,j,s,v} 22 6 P2(1,f)],sqrt(Pi(s,1))*(1-beta)*(B{s,i})',1);
            
            lmiterm([LMI{f,i,j,s,v} 21 7 P1(2,f)],sqrt(Pi(s,2))*beta*(B{s,i})',1);
            lmiterm([LMI{f,i,j,s,v} 21 7 X(v)],-sqrt(Pi(s,2))*beta,(B{s,i})');
            lmiterm([LMI{f,i,j,s,v} 21 8 P2(2,f)],sqrt(Pi(s,2))*beta*(B{s,i})',1);
            lmiterm([LMI{f,i,j,s,v} 22 7 P1(2,f)],sqrt(Pi(s,2))*(1-beta)*(B{s,i})',1);
            lmiterm([LMI{f,i,j,s,v} 22 7 X(v)],-sqrt(Pi(s,2))*(1-beta),(B{s,i})');
            lmiterm([LMI{f,i,j,s,v} 22 8 P2(2,f)],sqrt(Pi(s,2))*(1-beta)*(B{s,i})',1);
            
            lmiterm([LMI{f,i,j,s,v} 21 9 P1(3,f)],sqrt(Pi(s,3))*beta*(B{s,i})',1);
            lmiterm([LMI{f,i,j,s,v} 21 9 X(v)],-sqrt(Pi(s,3))*beta,(B{s,i})');
            lmiterm([LMI{f,i,j,s,v} 21 10 P2(3,f)],sqrt(Pi(s,3))*beta*(B{s,i})',1);
            lmiterm([LMI{f,i,j,s,v} 22 9 P1(3,f)],sqrt(Pi(s,3))*(1-beta)*(B{s,i})',1);
            lmiterm([LMI{f,i,j,s,v} 22 9 X(v)],-sqrt(Pi(s,3))*(1-beta),(B{s,i})');
            lmiterm([LMI{f,i,j,s,v} 22 10 P2(3,f)],sqrt(Pi(s,3))*(1-beta)*(B{s,i})',1);
            
            %omga66
            lmiterm([LMI{f,i,j,s,v} 21 21 X(v)],aib2,-1,'s');
            lmiterm([LMI{f,i,j,s,v} 22 22 X(v)],aib2,-1,'s');
           
    
 %ji%%%%%%%%%%%%%%%%%%%%%%%%%%
            lmiterm([LMI{f,i,j,s,v} 1 1 W1(s,v,j)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 1 2 W2(s,v,j)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 2 2 W3(s,v,j)],-1,1);
      
            lmiterm([LMI{f,i,j,s,v} 3 1 0],-S'*E{s,j});
            
            lmiterm([LMI{f,i,j,s,v} 3 3 0],-S'*F{s,j});
            lmiterm([LMI{f,i,j,s,v} 3 3 0],-F{s,j}'*S);
            lmiterm([LMI{f,i,j,s,v} 3 3 0],-R);
            lmiterm([LMI{f,i,j,s,v} 3 3 gamma],1,1);
            
            lmiterm([LMI{f,i,j,s,v} 4 1 0],Q1_*E{s,j});
            lmiterm([LMI{f,i,j,s,v} 4 3 0],Q1_*F{s,j});
            lmiterm([LMI{f,i,j,s,v} 4 4 0],-1);
            
            %omga41
            lmiterm([LMI{f,i,j,s,v} 5 1 P1(1,f)],1,sqrt(Pi(s,1))*A{s,j});
            lmiterm([LMI{f,i,j,s,v} 5 1 P2(1,f)],1,sqrt(Pi(s,1))*alpha*C{s,i});
            lmiterm([LMI{f,i,j,s,v} 5 2 P2(1,f)],1,sqrt(Pi(s,1))*(1-alpha));
            lmiterm([LMI{f,i,j,s,v} 6 1 -P2(1,f)],1,sqrt(Pi(s,1))*A{s,j});
            lmiterm([LMI{f,i,j,s,v} 6 1 P3(1,f)],1,sqrt(Pi(s,1))*alpha*C{s,i});
            lmiterm([LMI{f,i,j,s,v} 6 2 P3(1,f)],1,sqrt(Pi(s,1))*(1-alpha));           
            lmiterm([LMI{f,i,j,s,v} 5 1 Y(v)],sqrt(Pi(s,1))*beta*B{s,j},C{s,i});
            lmiterm([LMI{f,i,j,s,v} 5 2 Y(v)],sqrt(Pi(s,1))*(1-beta)*B{s,j},1);
            
            lmiterm([LMI{f,i,j,s,v} 7 1 P1(2,f)],1,sqrt(Pi(s,2))*A{s,j});
            lmiterm([LMI{f,i,j,s,v} 7 1 P2(2,f)],1,sqrt(Pi(s,2))*alpha*C{s,i});
            lmiterm([LMI{f,i,j,s,v} 7 2 P2(2,f)],1,sqrt(Pi(s,2))*(1-alpha));
            lmiterm([LMI{f,i,j,s,v} 8 1 -P2(2,f)],1,sqrt(Pi(s,2))*A{s,j});
            lmiterm([LMI{f,i,j,s,v} 8 1 P3(2,f)],1,sqrt(Pi(s,2))*alpha*C{s,i});
            lmiterm([LMI{f,i,j,s,v} 8 2 P3(2,f)],1,sqrt(Pi(s,2))*(1-alpha));
            lmiterm([LMI{f,i,j,s,v} 7 1 Y(v)],sqrt(Pi(s,2))*beta*B{s,j},C{s,i});
            lmiterm([LMI{f,i,j,s,v} 7 2 Y(v)],sqrt(Pi(s,2))*(1-beta)*B{s,j},1);
            
            lmiterm([LMI{f,i,j,s,v} 9 1 P1(3,f)],1,sqrt(Pi(s,3))*A{s,j});
            lmiterm([LMI{f,i,j,s,v} 9 1 P2(3,f)],1,sqrt(Pi(s,3))*alpha*C{s,i});
            lmiterm([LMI{f,i,j,s,v} 9 2 P2(3,f)],1,sqrt(Pi(s,3))*(1-alpha));
            lmiterm([LMI{f,i,j,s,v} 10 1 -P2(3,f)],1,sqrt(Pi(s,3))*A{s,j});
            lmiterm([LMI{f,i,j,s,v} 10 1 P3(3,f)],1,sqrt(Pi(s,3))*alpha*C{s,i});
            lmiterm([LMI{f,i,j,s,v} 10 2 P3(3,f)],1,sqrt(Pi(s,3))*(1-alpha));
            lmiterm([LMI{f,i,j,s,v} 9 1 Y(v)],sqrt(Pi(s,3))*beta*B{s,j},C{s,i});
            lmiterm([LMI{f,i,j,s,v} 9 2 Y(v)],sqrt(Pi(s,3))*(1-beta)*B{s,j},1);
            
            %omga42
            lmiterm([LMI{f,i,j,s,v} 5 3 P1(1,f)],1,sqrt(Pi(s,1))*D{s,j});
            lmiterm([LMI{f,i,j,s,v} 6 3 -P2(1,f)],1,sqrt(Pi(s,1))*D{s,j});
            
            lmiterm([LMI{f,i,j,s,v} 7 3 P1(2,f)],1,sqrt(Pi(s,2))*D{s,j});
            lmiterm([LMI{f,i,j,s,v} 8 3 -P2(2,f)],1,sqrt(Pi(s,2))*D{s,j});
            
            lmiterm([LMI{f,i,j,s,v} 9 3 P1(3,f)],1,sqrt(Pi(s,3))*D{s,j});
            lmiterm([LMI{f,i,j,s,v} 10 3 -P2(3,f)],1,sqrt(Pi(s,3))*D{s,j});
            
            %-omga44
            lmiterm([LMI{f,i,j,s,v} 5 5 P1(1,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 5 6 P2(1,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 6 6 P3(1,f)],-1,1);
            
            lmiterm([LMI{f,i,j,s,v} 7 7 P1(2,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 7 8 P2(2,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 8 8 P3(2,f)],-1,1);
            
            lmiterm([LMI{f,i,j,s,v} 9 9 P1(3,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 9 10 P2(3,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 10 10 P3(3,f)],-1,1);
            
            %omga51
            lmiterm([LMI{f,i,j,s,v} 11 1 Y(v)],sqrt(Pi(s,1))*beta_*B{s,j},C{s,i});
            lmiterm([LMI{f,i,j,s,v} 11 2 Y(v)],-sqrt(Pi(s,1))*beta_*B{s,j},1);
            
            
            lmiterm([LMI{f,i,j,s,v} 11 11 P1(1,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 11 12 P2(1,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 12 12 P3(1,f)],-1,1);
            
            lmiterm([LMI{f,i,j,s,v} 13 1 Y(v)],sqrt(Pi(s,2))*beta_*B{s,j},C{s,i});
            lmiterm([LMI{f,i,j,s,v} 13 2 Y(v)],-sqrt(Pi(s,2))*beta_*B{s,j},1);
            
            lmiterm([LMI{f,i,j,s,v} 13 13 P1(2,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 13 14 P2(2,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 14 14 P3(2,f)],-1,1);
            
            lmiterm([LMI{f,i,j,s,v} 15 1 Y(v)],sqrt(Pi(s,3))*beta_*B{s,j},C{s,i});
            lmiterm([LMI{f,i,j,s,v} 15 2 Y(v)],-sqrt(Pi(s,3))*beta_*B{s,j},1);
            
            lmiterm([LMI{f,i,j,s,v} 15 15 P1(3,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 15 16 P2(3,f)],-1,1);
            lmiterm([LMI{f,i,j,s,v} 16 16 P3(3,f)],-1,1);
            
            lmiterm([LMI{f,i,j,s,v} 17 1 0],kesi3{v}*C{s,i});
            lmiterm([LMI{f,i,j,s,v} 17 17 0],-aib1);
            lmiterm([LMI{f,i,j,s,v} 18 2 0],kesi3{v});
            lmiterm([LMI{f,i,j,s,v} 18 18 0],-aib1);
            
            lmiterm([LMI{f,i,j,s,v} 19 5 P1(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*beta*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 19 6 P2(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*beta*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 20 5 P1(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*(1-beta)*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 20 6 P2(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*(1-beta)*B{s,j}',1);
            
            lmiterm([LMI{f,i,j,s,v} 19 7 P1(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*beta*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 19 8 P2(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*beta*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 20 7 P1(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*(1-beta)*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 20 8 P2(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*(1-beta)*B{s,j}',1);
            
            lmiterm([LMI{f,i,j,s,v} 19 9 P1(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*beta*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 19 10 P2(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*beta*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 20 9 P1(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*(1-beta)*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 20 10 P2(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*(1-beta)*B{s,j}',1);
            
            lmiterm([LMI{f,i,j,s,v} 19 11 P1(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*beta_*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 19 12 P2(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*beta_*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 20 11 P1(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*(-beta_)*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 20 12 P2(1,f)],aib1*sqrt(Pi(s,1))*kesi1{v}*(-beta_)*B{s,j}',1);
            
            lmiterm([LMI{f,i,j,s,v} 19 13 P1(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*beta_*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 19 14 P2(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*beta_*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 20 13 P1(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*(-beta_)*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 20 14 P2(2,f)],aib1*sqrt(Pi(s,2))*kesi1{v}*(-beta_)*B{s,j}',1);
            
            lmiterm([LMI{f,i,j,s,v} 19 15 P1(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*beta_*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 19 16 P2(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*beta_*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 20 15 P1(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*(-beta_)*B{s,j}',1);
            lmiterm([LMI{f,i,j,s,v} 20 16 P2(3,f)],aib1*sqrt(Pi(s,3))*kesi1{v}*(-beta_)*B{s,j}',1);
            
            lmiterm([LMI{f,i,j,s,v} 19 19 0],-aib1);
            lmiterm([LMI{f,i,j,s,v} 20 20 0],-aib1);
            
            
            lmiterm([LMI{f,i,j,s,v} 21 1 Y(v)],aib2,C{s,i});
            lmiterm([LMI{f,i,j,s,v} 22 2 Y(v)],aib2,1);
            
            lmiterm([LMI{f,i,j,s,v} 21 5 P1(1,f)],sqrt(Pi(s,1))*beta*(B{s,j})',1);
            lmiterm([LMI{f,i,j,s,v} 21 5 X(v)],-sqrt(Pi(s,1))*beta,(B{s,j})');
            lmiterm([LMI{f,i,j,s,v} 21 6 P2(1,f)],sqrt(Pi(s,1))*beta*(B{s,j})',1);
            lmiterm([LMI{f,i,j,s,v} 22 5 P1(1,f)],sqrt(Pi(s,1))*(1-beta)*(B{s,j})',1);
            lmiterm([LMI{f,i,j,s,v} 22 5 X(v)],-sqrt(Pi(s,1))*(1-beta),(B{s,j})');
            lmiterm([LMI{f,i,j,s,v} 22 6 P2(1,f)],sqrt(Pi(s,1))*(1-beta)*(B{s,j})',1);
            
            lmiterm([LMI{f,i,j,s,v} 21 7 P1(2,f)],sqrt(Pi(s,2))*beta*(B{s,j})',1);
            lmiterm([LMI{f,i,j,s,v} 21 7 X(v)],-sqrt(Pi(s,2))*beta,(B{s,j})');
            lmiterm([LMI{f,i,j,s,v} 21 8 P2(2,f)],sqrt(Pi(s,2))*beta*(B{s,j})',1);
            lmiterm([LMI{f,i,j,s,v} 22 7 P1(2,f)],sqrt(Pi(s,2))*(1-beta)*(B{s,j})',1);
            lmiterm([LMI{f,i,j,s,v} 22 7 X(v)],-sqrt(Pi(s,2))*(1-beta),(B{s,j})');
            lmiterm([LMI{f,i,j,s,v} 22 8 P2(2,f)],sqrt(Pi(s,2))*(1-beta)*(B{s,j})',1);
            
            lmiterm([LMI{f,i,j,s,v} 21 9 P1(3,f)],sqrt(Pi(s,3))*beta*(B{s,j})',1);
            lmiterm([LMI{f,i,j,s,v} 21 9 X(v)],-sqrt(Pi(s,3))*beta,(B{s,j})');
            lmiterm([LMI{f,i,j,s,v} 21 10 P2(3,f)],sqrt(Pi(s,3))*beta*(B{s,j})',1);
            lmiterm([LMI{f,i,j,s,v} 22 9 P1(3,f)],sqrt(Pi(s,3))*(1-beta)*(B{s,j})',1);
            lmiterm([LMI{f,i,j,s,v} 22 9 X(v)],-sqrt(Pi(s,3))*(1-beta),(B{s,j})');
            lmiterm([LMI{f,i,j,s,v} 22 10 P2(3,f)],sqrt(Pi(s,3))*(1-beta)*(B{s,j})',1);
            
            %omga66
            lmiterm([LMI{f,i,j,s,v} 21 21 X(v)],aib2,-1,'s');
            lmiterm([LMI{f,i,j,s,v} 22 22 X(v)],aib2,-1,'s');
            
                end 
            end
        end
    end
end

%% solve LMIs
lmisys=getlmis;
y=decnbr(lmisys);              %得到lmis系统变量个数
c=[zeros(y-1,1);-1]; 
options=[1e-5,1000,0,0,0];     % 精度要求为1e-5，最大迭代次数为1000
[a1,b1]=mincx(lmisys,c,options);
gamma=dec2mat(lmisys,b1,gamma)

for v=1:S2                                %控制器
    
           fprintf('K_%d%d  is',v)              % %d%d两个输出之间没有逗号间隔
          
          dec2mat(lmisys,b1,X(v))^(-1)*(dec2mat(lmisys,b1,Y(v)))        %dec2mat：决策变量转化为矩阵形式 
end       

toc      %计时



