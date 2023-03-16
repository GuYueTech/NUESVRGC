%% benchmark network given
 N=30; % number of nodes in network
avg_degree=3; 


% A= random_net(N, avg_degree);
% A= scalefree_BA(N, avg_degree);
%
p=0.2; A= smallworld_net(N, avg_degree, p);
% A=A;

% load SW100conn.mat;
% A=B;

w=-2 + (4).*rand(N,1);
n_time_points=2;
n_expers=4000;
noise=0.01;
c=1;

XX = kuramoto_dis(w,A,c,n_time_points,n_expers,noise);



% [1,2,3,4,5]*sin([2;3;4;5;6]-1)