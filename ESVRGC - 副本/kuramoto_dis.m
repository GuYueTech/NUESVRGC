
% % %% benchmark network given
% %  N=10; % number of nodes in network
% % avg_degree=3; 
% % 
% % 
% %  %A= random_net(N, avg_degree);
% % % A= scalefree_BA(N, avg_degree);
% % p=0.2; A= smallworld_net(N, avg_degree, p);
% % 
% % 
% % % load SW100conn.mat;
% % % A=B;
% % 
% % w=-2 + (4).*rand(N,1);
% % n_time_points=
% % n_expers=
% % noise=

function XX= kuramoto_dis(w,A,c,n_time_points,n_expers,noise)

N=size(A,1);
T=30;
X=zeros(N,T); %records of states: T columns
% rng(1,'v5uniform')
% X(:,1)=rand(N,1); %set initial state
samplerate = 1;
% n_expers=500;
%---Parameters for gene regulatory dynamic-----------


XX=[];
for m=1:n_expers
%  X(:,1)=log(m)*rand(N,1); %set initial state   
% X(:,1)=rand(N,1); %set initial state   
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
X(:,1)=sqrt(m)*rand(N,1); %set initial state  
%X(:,1)=-3.14 +(3.14+3.14)*rand(N,1)+1/m; 
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% R=noise*randn(N,1);
for i=1:N
for k=2:T %real time from 0:1:50
     X(i,k)=X(i,k-1) +samplerate*w(i) +samplerate*c*A(i,:)*sin(X(:,k-1)-X(i,k-1)) + noise*randn;
%     X(i,k)=X(i,k-1) +w(i) +samplerate*A(i,:)*sin(X(:,k-1)-X(i,k-1)) + noise*randn;
end	  
end
XT=X(:,1:n_time_points);
   XX=[XX XT];
   

end


   figure; 

plot(X'); 







