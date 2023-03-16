
clc;
clear all;
close all;


%% 因果性推断 网络重构 %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
    

M=6; %number of simulated time series
numsimu=1; % number of simulations 
alpha=0.01; % P-value for GC F-test

	
B= [ 0 0 0 0 1 0 ;
     0 0 0 1 0 1 ;
     0 0 0 0 1 0 ;
     1 0 1 0 0 0 ;
     0 1 0 0 0 0 ;
     1 0 1 0 0 0];	


%% -----------------generate data ---------------------------



FSNCGC_AUROC=zeros(numsimu,1);
FSNCGC_AUPR=zeros(numsimu,1);




  
 for i=1:numsimu
    disp(['Simulation ' int2str(i) ' of ' int2str(numsimu)]);
%% Generating Data    
 
% %% Generating Data    
%      n_timepoints = 30; %时间点要比maxlag大
%      m_delta=30;
%     noise=0.01; maxlag=10;
%     % x = MG_dis_my(B,noise,n_timepoints,n_measurements);
%     x = MG_dis_my_maxlag(B,noise,n_timepoints,m_delta,maxlag);
%     Y=x';
   


t1=cputime;

    n_timepoints = 30; %时间点要比maxlag大，还有maxorder大
    m_delta=10;
    [XX]=GRN_con_my_maxlag_final(B,n_timepoints, m_delta);
    Y=XX;
%     Yp=XXp;
    

fprintf('calculate the marginal dmi costs %5.1fs.\n', cputime-t1);
  

    
     
     %% KCV for model order selection
%     kfold=5; %must greater than one;
%     maxorder=20; %不能大于 n_timepoints
%     z1=estimate_order_CGC( Y,kfold, maxorder,n_timepoints);
%     pmax=z1.order; 
    pmax=10;
     
     
%%  estimation of NCGC
 %---------------rbf center selection ----------cross validation-----------------   
%    kfold=5; %must greater than one;
%    maxCenter=10;
%    [z, C]= estimate_numC_NCGC( Y, kfold, maxCenter,pmax ,n_timepoints);

%     rbf_C=z.C; 
%     num_C=size(rbf_C,1);   
    
      
      
      
      
%%  estimation of FSNCGC
% 

    rbf_C=[]; 
    num_C=5;  
      [NENCGC1,pNENCGC1,NErbf_Center32,NEind3,NEPhi_g32,NEecC1,NEecC2]=NENCGC(Y,n_timepoints,pmax,num_C,rbf_C,1);%数据标准化了      
      NENCGC1= NENCGC1.*(ones(M,M)-eye(M,M));
      pNENCGC1=pNENCGC1.*(ones(M,M)-eye(M,M));  
      
      
    
      FSNCGC1= NENCGC1;
      pFSNCGC1=pNENCGC1;

      
      
      
      
      
      

%% COMPUTE AUROC & AUPR
%不设置pval截断，直接用推得的矩阵
% CGC1=CGC1.*(pCGC1<alpha);
% GLasso_CGC1=GLasso_CGC1.*(pGLasso_CGC1<alpha);
% mBTSCGC1=mBTSCGC1.*(pmBTSCGC1<alpha);
% NCGC1=NCGC1.*(pNCGC1<alpha);
% NENCGC1=NCGC1.*(pNENCGC1<alpha);



%FSNCGC
FSNCGC1=(FSNCGC1-min(min(FSNCGC1)))/(max(max(FSNCGC1))-min(min(FSNCGC1)));
[FSNCGC_AUROC1, ~ ,FSNCGC_AUPR1]= AUCs(FSNCGC1, B);


FSNCGC_AUROC(i)=FSNCGC_AUROC1;
FSNCGC_AUPR(i)=FSNCGC_AUPR1;

FSNCGCM(:,:,i)=FSNCGC1; pFSNCGCM(:,:,i)=pFSNCGC1;  

%% COMPUTE other metrics
%pval截断，直接用推得的矩阵

FSNCGC2=FSNCGC1.*(pFSNCGC1<alpha);
[~,~,~,~,FSNCGC_TPR(i), FSNCGC_FPR(i), FSNCGC_PPV(i), FSNCGC_Fscore(i), FSNCGC_ACC(i), FSNCGC_MCC(i)] = cal_statistics(FSNCGC2, B);

end




Aver_FSNCGC_AUROC=sum(FSNCGC_AUROC)/numsimu;      Std_FSNCGC_AUROC=std(FSNCGC_AUROC);
Aver_FSNCGC_AUPR=sum(FSNCGC_AUPR)/numsimu;        Std_FSNCGC_AUPR=std(FSNCGC_AUPR);
Aver_FSNCGC_TPR=sum(FSNCGC_TPR)/numsimu;          Std_FSNCGC_TPR=std(FSNCGC_TPR);
Aver_FSNCGC_FPR=sum(FSNCGC_FPR)/numsimu;          Std_FSNCGC_FPR=std(FSNCGC_FPR);
Aver_FSNCGC_PPV=sum(FSNCGC_PPV)/numsimu;          Std_FSNCGC_PPV=std(FSNCGC_PPV);
Aver_FSNCGC_Fscore=sum(FSNCGC_Fscore)/numsimu;    Std_FSNCGC_Fscore=std(FSNCGC_Fscore);
Aver_FSNCGC_ACC=sum(FSNCGC_ACC)/numsimu;          Std_FSNCGC_ACC=std(FSNCGC_ACC);
Aver_FSNCGC_MCC=sum(FSNCGC_MCC)/numsimu;          Std_FSNCGC_MCC=std(FSNCGC_MCC);






%%		
Result1=[
         Aver_FSNCGC_AUROC         Std_FSNCGC_AUROC; 
         Aver_FSNCGC_AUPR          Std_FSNCGC_AUPR; 
         Aver_FSNCGC_TPR           Std_FSNCGC_TPR; 
         Aver_FSNCGC_FPR           Std_FSNCGC_FPR; 
         Aver_FSNCGC_PPV           Std_FSNCGC_PPV; 
%          Aver_FSNCGC_Fscore        Std_FSNCGC_Fscore; 
         Aver_FSNCGC_ACC           Std_FSNCGC_ACC;
%          Aver_FSNCGC_MCC           Std_FSNCGC_MCC
]; 
%%		 
Result=Result1';		 
		 
ind=0;
for i=1:M
    for j=1:M
        if i~=j
    ind=ind+1;
    pFSNCGC_sum(i,j)=sum(pFSNCGCM(i,j,:)<alpha);
        end
    end
end



%% plot true matrix
figure;subplot(121);
imagesc(B); 
set(gca,'XTickLabel',{}); 
set(gca,'YTickLabel',{});
title('Connectivity','FontName','Times New Roman','FontSize',14);

%% plot number of discovery rate matrix
subplot(122);
imagesc(pFSNCGC_sum);
% set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10])% X坐标轴刻度数据点位置  
set(gca,'XTickLabel',{}) %X坐标轴刻度处显示的字符 
set(gca,'FontName','Times New Roman','FontSize',14)%设置坐标轴刻度字体名

% set(gca, 'YTick', [1 2 3 4 5 6 7 8 9 10])% X坐标轴刻度数据点位置  
set(gca,'YTickLabel',{}) %X坐标轴刻度处显示的字符 
set(gca,'FontName','Times New Roman','FontSize',14)%设置坐标轴刻度字体名
title('FSNCGC','FontName','Times New Roman','FontSize',14);    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    





%% %% 步长是变的，dde23是变步长的solver
% figure;
% plot(sol.x,sol.y(1,:),'-r' )
% hold on
% plot(sol.x,sol.y(2,:),'-b' )
% hold on;
% plot(sol.x,sol.y(3,:),'-g' )
% hold off;
% title('Mackey-Glass');
% xlabel('t');
% ylabel('x');
% legend('x_1','x_2','x_3');
% 
% figure;
% plot(sol.x,sol.yp(1,:),'r-o' )
% hold on
% plot(sol.x,sol.yp(2,:),'b-o' )
% hold on;
% plot(sol.x,sol.yp(3,:),'g-o' )
% hold off;
% title('Mackey-Glass');
% xlabel('t');
% ylabel('${\dot x}$','Interpreter','latex');
% legend({'${\dot x}_1$','${\dot x}_2$','${\dot x}_3$'},'Interpreter','latex');
% %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
% 
% 
% 
% %[data_new, derivative, index] = estimatediff(data, t, type, kth, order);
%  k_order=1;
% [data_new, derivative]=estimatediff(sol.y', sol.x', 'solver', k_order, []);
% which_var=1;
% derivative_true=sol.yp';
%     figure; 
%     plot(data_new(:,which_var),'r');
%     figure; 
%     plot(derivative_true(:,which_var),'r-o');
%     hold on; 
%     plot(derivative(:,which_var),'b-*'); 
%%