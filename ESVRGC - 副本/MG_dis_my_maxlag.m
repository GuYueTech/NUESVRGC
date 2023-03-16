


function [XX] = MG_dis_my_maxlag(A,noise,n_timepoints,M,maxlag)
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



N=size(A,1);
T=300;
X=zeros(N,T); %records of states: T rows
% rng(1,'v5uniform')
% X(:,1)=rand(N,1); %set initial state
samplerate = 1;
% M=3;

XX=[];
for m=1:M
%  X(:,1)=log(m)*rand(N,1); %set initial state   
% X(:,1)=rand(N,1); %set initial state   
X(:,1:maxlag)=sqrt(m)*rand(N,maxlag); %set initial state   

% X(:,1:10)=ones(N,10); %set initial state   
for k=maxlag+1:T %real time from 0:1:50
    X(:,k)=X(:,k-1) -samplerate*0.1*X(:,k-1)  +samplerate*A*(X(:,k-maxlag)./(X(:,k-maxlag).^10+1)) + noise*randn(N,1);
   
end
XT=X(:,maxlag+1:maxlag+n_timepoints);
   XX=[XX XT];
end

% figure; 
% 
% plot(X'); 












