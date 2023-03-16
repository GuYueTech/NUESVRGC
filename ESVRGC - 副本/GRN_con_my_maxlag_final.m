


%function [XX, XXp] = MG_con_my_maxlag(A,n_timepoints,M)
%% 时滞为 1
%%
% A---邻接矩阵
% noise---噪声
% n_timepoints----每次测量的时间点个数
% M-------测量次数
%% benchmark network given
%  N=100; % number of nodes in network
%  avg_degree=5; 
%  A= random_net(N, avg_degree);
%  A= scalefree_BA(N, avg_degree);
%  p=0.2; A= smallworld_net(N, avg_degree, p);
%  load SW100conn.mat;
%  A=B;


% noise=0.01;
% A=[ 0,0,0,0,1
%     1,0,1,0,0
%     1,0,0,0,0
%     1,0,0,0,0
%     0,0,0,1,0];
% % A=[ 0,1,0,0,0,0;
% %     0,1,1,0,0,0;
% %     0,1,0,0,0,1;
% %     0,0,0,1,1,1;
% %     0,0,0,1,0,0;
% %     0,0,0,0,1,0];



% N=size(A,1);
% T=100;
% X=zeros(N,T); %records of states: T columns
% % rng(1,'v5uniform')
% % X(:,1)=rand(N,1); %set initial state
% samplerate = 1;
% % M=3;
% 
% XX=[];
% for m=1:M
% 
% X(:,1:maxlag)=sqrt(m)*rand(N,maxlag); %set initial state   
% 
% % X(:,1:10)=ones(N,10); %set initial state   
% for k=maxlag+1:T %real time from 0:1:50
%     X(:,k)=X(:,k-1) -samplerate*0.1*X(:,k-1)  +samplerate*A*(X(:,k-maxlag)./(X(:,k-maxlag).^10+1)) + noise*randn(N,1);
%    
% end
% XT=X(:,maxlag+1:maxlag+n_timepoints);
%    XX=[XX XT];
% end



% function [XX, XXp] = GRN_con_my_maxlag_final(A,n_timepoints,M)
function [XX] = GRN_con_my_maxlag_final(A,n_timepoints,M)
%% 时滞为 1
%%
% A---邻接矩阵
% noise---噪声
% n_timepoints----每次测量的时间点个数
% M-------测量次数

%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%


% clc
% clear
% 
% A=[ 0 0 1;
%     1 0 0;
%     0 1 0];
% n_timepoints=30;
% M=30;



%% 生成数据 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
N=size(A,1);
noise=0.01;
n=1;
beta=1;
XX=[];
XXp=[];

for m=1:M
history = @(t) (m)*rand(N,1);
% history = @(t) sqrt(m)*rand(N,1);
ddefun = @(t,y,Z) [ -beta*y(1)+Z(5,1)^n/(Z(5,1)^n+1);
                    -beta*y(2)+Z(4,1)^n/(Z(4,1)^n+1)+Z(6,1)^n/(Z(6,1)^n+1);
                    -beta*y(3)+Z(5,1)^n/(Z(5,1)^n+1);
					-beta*y(4)+Z(1,1)^n/(Z(1,1)^n+1)+Z(3,1)^n/(Z(3,1)^n+1);
					-beta*y(5)+Z(2,1)^n/(Z(2,1)^n+1);
                    -beta*y(6)+Z(1,1)^n/(Z(1,1)^n+1)+Z(3,1)^n/(Z(3,1)^n+1)]+ noise*randn(N,1);
lags=6;
%tspan=[0,100];
t0=0;
tf=59;
tspan=[t0,tf]; 
sol = dde23(ddefun,lags,history,tspan);

%% 等间隔步长取值 

dt=0.01;
tint1 = t0:dt:tf;
[yint1,yintp1] = deval(sol,tint1);


%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
XT=yint1(:,(lags+1)/dt:1/dt:((lags+1)/dt+n_timepoints/dt));
XX=[XX XT];

%%假设导数可测量到，得出如下式子
% XTp=yintp1(:,(lags+1)/dt:1/dt:((lags+1)/dt+n_timepoints/dt));
% XXp=[XXp XTp];
%% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%




%% 因为时延 tau, 所以数据取 tau 以后的值
%% for 导数估计
tint=tint1((lags+1)/dt:end);
yint=yint1(:,(lags+1)/dt:end);
yintp=yintp1(:,(lags+1)/dt:end);

% subplot(121);
% plot(tint,yint);
% subplot(122);
% plot(tint,yintp);


%% 数据曲线，各个采样点的数据
% figure;
% %plot(tint,yint);
% plot(tint,yint(1,:),'-m' );
% hold on
% plot(tint,yint(2,:),'-b' );
% hold on;
% plot(tint,yint(3,:),'-g' );
% hold on;
% plot(tint,yint(4,:),'-r' );
% hold on;
% plot(tint,yint(5,:),'-y' );
% hold on;
% plot(tint,yint(6,:),'-c' );
% hold off;
% title('GRN');
% xlabel('t','Interpreter','latex'); 
% ylabel('$x$','Interpreter','latex');
% legend({'$x_1$','$x_2$','$x_3$','$x_4$','$x_5$','$x_6$'},'Interpreter','latex');% 注意和matlab中解释器‘tex’的区别 是需要加上$$么？




% 
% %% 导数不可测量，先用估计的方法去求
% [data_new, derivative,  derivative_true, index] = estimatediff_my_final(yint', tint',yintp', 'solver',1, []);
% % [data_new, derivative,  derivative_true, index] = estimatediff_my(yint', tint',yintp','eularforward');
% 
% 
% % XT=data_new(lags/dt: 1/(dt):(lags/dt+(n_timepoints-1)/(dt)),:);  %1/(dt)
% % XX=[XX; XT];
% % 
% % XTp=derivative(lags/dt: 1/(dt):(lags/dt+(n_timepoints-1)/(dt)),:);
% % XXp=[XXp; XTp];
% 
% XT=data_new((lags+1)/dt: 1/(dt):((lags+1)/dt+(n_timepoints-1)/(dt)),:);  %1/(dt)
% XX=[XX; XT];
% 
% XTp=derivative((lags+1)/dt: 1/(dt):((lags+1)/dt+(n_timepoints-1)/(dt)),:);
% XXp=[XXp; XTp];
% 
% % XT=data_new(1: 1/(dt):(1+n_timepoints/(dt)),:);  %1/(dt)
% % XX=[XX; XT];
% % 
% % XTp=derivative(1: 1/(dt):(1+n_timepoints/(dt)),:);
% % XXp=[XXp; XTp];


end








