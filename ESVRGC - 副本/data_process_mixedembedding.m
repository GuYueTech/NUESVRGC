function [y_S,Phi_g2,Phi_g3,ind3,C,ecC1,ecC2]=data_process_mixedembedding(data,n_time_points,D,num_C,rbf_C,nonlinear)
%% 非线性加入


%%%%% input %%%%%%%%%%%
%data is n_samples*nvars 样本*变量数
%n_time_points in each experiment
%D-1 is model order
%nonlinear is true ---> kernel matrix of var  
%%%%% output %%%%%%%%%%%%
% y_S : target variable matrix m_Delta*(T-D+1) * nvar ; rows:m_Delta个(D:end)行
% Phi_g2 : matrix of var
% Phi_g3 : matrix of nonlinear var
% C : Centers by Kmeans or Fuzzy C-Means


[n_samples, N]=size(data);
%1
% av = sum(data,1)/n_samples;
% data=data - ones(size(data)) * diag(av);
% % Phi_g1=Phi_g1./repmat(std(Phi_g1),size(Phi_g1,1),1);
% data= normc(data); 

%2
% minallV1=min(data);
% rang1=kron((1./range(data)),ones(n_samples,1));
% data=(data-kron(minallV1,ones(n_samples,1))).*rang1;

if size(n_time_points,2)~=1
%% 时间点在每次实验当中都不一样，此时n_time_points是一个矢量
m_Delta=size(n_time_points,2);
%%每次试验的数据分别存入相应的cell中  
tt=0;
for i = 1:m_Delta
         temp = data((tt+1):(tt+n_time_points(i)),:);
         data_cell{1,i} = temp;  
		 tt=tt+n_time_points(i);
end
S_target=1:N; 
%%%%%%%% Forming y and Phi and Preprocessing %%%%%%%%%%%
y_S = [];  
Phi_g1 = [];  
for j = 1:m_Delta
    Phi_g1_temp = [];
    y_S = [y_S;data_cell{1,j}(D:end,S_target)];
    for d = 1:D-1
        Phi_g1_temp = [Phi_g1_temp,data_cell{1,j}(D-d:end-d,:)];
    end
    Phi_g1 = [Phi_g1;Phi_g1_temp];
end

else


m_Delta=n_samples/n_time_points;% number of experiments/measurements
%%每次试验的数据分别存入相应的cell中  
for i = 1:m_Delta
         temp = data(((i-1)*n_time_points+1):(i*n_time_points),:);
         data_cell{1,i} = temp;  
end
S_target=1:N; 
%%%%%%%% Forming y and Phi and Preprocessing %%%%%%%%%%%
y_S = [];  
Phi_g1 = [];  
for j = 1:m_Delta
    Phi_g1_temp = [];
    y_S = [y_S;data_cell{1,j}(D:end,S_target)];
    for d = 1:D-1
        Phi_g1_temp = [Phi_g1_temp,data_cell{1,j}(D-d:end-d,:)];
    end
    Phi_g1 = [Phi_g1;Phi_g1_temp];
end
%%%% Phi_g1=[X1(D-1) X2(D-1)...XN(D-1)|...|X1(1) X2(1)...XN(1)]
end





%% %%%%%%%%%% Normalization %%%%%%%%%%%%%%%%%
[m,N_t] = size(Phi_g1);
% num_C=D-1; %number of centers

temp =[];
sigma=1;

av = sum(Phi_g1,1)/m;
Phi_g1 = Phi_g1 - ones(size(Phi_g1)) * diag(av);
% Phi_g1=Phi_g1./repmat(std(Phi_g1),size(Phi_g1,1),1);
Phi_g1 = normc(Phi_g1);  % We normalize columns of Phi

av2 = sum(y_S,1)/m;
y_S = y_S - ones(size(y_S)) * diag(av2);
% y_S=y_S./repmat(std(y_S),size(y_S,1),1);
y_S = normc(y_S);  % We normalize columns of Phi
%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%

% Standardization of the input matrix columnwise in [0,1].
% minallV=min(Phi_g1);
% rang=kron((1./range(Phi_g1)),ones(m,1));
% Phi_g1=(Phi_g1-kron(minallV,ones(m,1))).*rang;
% 
% minallV1=min(y_S);
% rang1=kron((1./range(y_S)),ones(m,1));
% y_S=(y_S-kron(minallV1,ones(m,1))).*rang1;



%% %%%%% rearrange the Phi_g1  %%%%%%%%%%%%
temp=[];
for i=1:N
    pool_Si = S_target(i)+(0:N:(D-2)*N);
    temp=[temp Phi_g1(:,pool_Si)];
end
    Phi_g2=temp;
%%%% Phi_g2=[X1(D-1) X1(D-2)...X1(1)|...|XN(D-1) XN(D-2)...XN(1)]
%%%% using Phi_g2 for kernel regression 
if nonlinear
%% mixed embedding   
    pmax=D-1;
   ecC1 = MixedEmbeddingMulti(data,pmax,1,5,0.02,2,n_time_points);   
% 将矩阵的行 按照第一列元素的大小 做由小到大排列 
  ecC2=cell(N,1);
  for KK=1:N
    [~,idx1]=sort(ecC1{KK,1}(1:end-1,1));
    ecC2{KK}=ecC1{KK,1}(idx1,:);
  end
%%    

  sigma=1;
  if isempty(rbf_C)
  [idx, C]=kmeans(Phi_g2(:,:),num_C,'emptyaction','drop');
% [C, ~,~] = FCMClust(Phi_g2,num_C);
    else
      C=rbf_C;
  end
 
  
 Phi_g3=cell(N,1); ind3=cell(N,1);
  for  j=1:N
  temp2=[];
   for i=1:m
    temp1=[];  indk=[];
    for k=1:N
	     idx2=find(ecC2{j}(:,1)==k);
        if ~isempty(idx2)
          klags=ecC2{j}(idx2,2);
		  indk=[indk k];
            for i_C=1:num_C
%% consider only mixed embedding lags as inputs


 %         temp1=[temp1
 %         exp(-sum((Phi_g2(i,klags)-C(i_C,klags)).^2)/(2*sigma*sigma))]; %
 %         有问题 叉叉。。。。   之前出现的低秩矩阵的问题 就是这里弄错了, 改成下面的就好了哦
        temp1=[temp1 exp(-sum((Phi_g2(i,(k-1)*(D-1)+klags)-C(i_C,(k-1)*(D-1)+klags)).^2)/(2*sigma*sigma))]; 
 
%         temp1=[temp1 exp(-sum(((Phi_g2(i,((k-1)*(D-1)+1):(k*(D-1)))-C(i_C,((k-1)*(D-1)+1):(k*(D-1))))).^2)/(2*sigma*sigma))];      
            end
        end
    end
       temp2=[temp2;temp1];
    end
    Phi_g3{j}=temp2; 
	ind3{j}=indk;
   av3 = sum(Phi_g3{j},1)/ size(Phi_g3{j},1);
   Phi_g3{j} = Phi_g3{j} - ones(size(Phi_g3{j})) * diag(av3);
% % Phi_g1=Phi_g1./repmat(std(Phi_g1),size(Phi_g1,1),1);
   Phi_g3{j} = normc(Phi_g3{j});  % We normalize columns of Phi
   end

 else
    Phi_g3=[];
    C=[];
end

%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
end