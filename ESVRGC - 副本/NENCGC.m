
function [GC,p_value,C,ind3,Phi_g3,ecC1,ecC2]=NENCGC(data,n_time_points,pmax,num_C,rbf_C,nonlinear)

% data： 行是样本，列是变量
% n_time_points = 20; 
% D=2;%time-lags, D-1 is the model order
% nonlinear=false;
% num_C=[];


% C : Centers by Kmeans or Fuzzy C-Means
D=pmax+1;
[y_S,~,Phi_g3,ind3,C,ecC1,ecC2]=data_process_mixedembedding(data,n_time_points,D,num_C,rbf_C,nonlinear);
%按照VAR的格式构造数据矩阵

N=size(data,2);
p_value=ones(N,N);
GC=zeros(N,N);

%[m,N_t] = size(Phi_g2);% 求出构造好的矩阵的行和列

for i=1:N
    Y=y_S(:,i);
	Phi_U=Phi_g3{i};

% Unrestricted regressions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    [SigmaU, Upu]=Regression_LS(Y,Phi_U);

    Nu=size(ind3{i},2)*num_C;
    
% Restricted regressions %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%  
 for ii=1:size(ind3{i},2)
     Phi_R=[];  
     Phi_R=Phi_g3{i};
     Phi_R(:,((ii-1)*num_C+1):ii*num_C)=[];
     [SigmaR, Upr]=Regression_LS(Y,Phi_R);
    
  GC(i,ind3{i}(ii))=log(SigmaR./SigmaU); 
  Nr=(size(ind3{i},2)-1)*num_C;
  
   p_value(i,ind3{i}(ii)) = egc_LinReg_Ftest(Upu',Upr',Nu,Nr);
 end
 end
% M = size(GC,2);
% GC = GC.*(ones(M,M)-eye(M,M));
% GC=(GC-min(min(GC)))/(max(max(GC))-min(min(GC)));

for i =1:1:size(GC,2)
    GC(i,i) = 0;
end
 end


 function [Sigma, Up]= Regression_LS(Y,Phi) % Linear Regression
 %% 
 %sigma : 误差方差
 %Up    ：回归误差
 %%
  if isempty(Phi) %if no conditioning, ce will be the Entropy of B
    Sigma=var(Y);
    Up=Y';
  else
     %%
    Yb=Y'; % inversion works with data organized in rows
    Z=Phi';
    Am=Yb/Z; % least squares!

    Yp=Am*Z; 
    Up=Yb-Yp;
    Sigma=cov(Up');
  %%
%       Z=inv(Phi'*Phi)*Phi'*Y;
%       Yp=Phi*Z;
%       Up=Y-Yp;
%       Sigma=cov(Up);

%%

%     Am=Phi\Y; % least squares!
% 
%     Yp=Phi*Am; 
%     Up=Y-Yp;
%     Sigma=cov(Up);  
  
    
  end

 end
 













