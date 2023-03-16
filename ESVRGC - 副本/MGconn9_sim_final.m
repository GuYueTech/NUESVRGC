clc;
clear all; 
close all; 
M=9; %number of simulated time series
numsimu=10; % number of simulations 
alpha=0.01; % P-value for GC F-test




	
B=[ 0 0 0 0 1 0 0 0 0 ;
    0 0 0 1 0 1 0 0 0 ;
    0 0 0 0 0 0 1 0 0 ;
    0 1 0 0 0 0 0 1 0 ;
    1 0 0 0 0 0 0 0 1 ;
    0 1 0 0 0 0 0 1 0 ;
    0 0 1 0 0 0 0 0 0 ;
    0 0 0 1 0 1 0 0 0 ;
    0 0 0 0 1 0 0 0 0 ];

% figure;
% imagesc(B); 
% set(gca,'XTickLabel',{}); 
% set(gca,'YTickLabel',{});
% title('Connectivity','FontName','Times New Roman','FontSize',14);





%% -----------------generate data ---------------------------


CGC_AUROC=zeros(numsimu,1);
CGC_AUPR=zeros(numsimu,1);
GLasso_CGC_AUROC=zeros(numsimu,1);
GLasso_CGC_AUPR=zeros(numsimu,1);
mBTSCGC_AUROC=zeros(numsimu,1);
mBTSCGC_AUPR=zeros(numsimu,1);
NCGC_AUROC=zeros(numsimu,1);
NCGC_AUPR=zeros(numsimu,1);
NENCGC_AUROC=zeros(numsimu,1);
NENCGC_AUPR=zeros(numsimu,1);


  
 for i=1:numsimu
    disp(['Simulation ' int2str(i) ' of ' int2str(numsimu)]);
%% Generating Data    
 
%% Generating Data    
     n_timepoints = 30; %时间点要比maxlag大
     m_delta=30;
    noise=0.01; maxlag=10;
    % x = MG_dis_my(B,noise,n_timepoints,n_measurements);
    x = MG_dis_my_maxlag(B,noise,n_timepoints,m_delta,maxlag);
     Y=x';
     

     
     %% KCV for model order selection
    kfold=5; %must greater than one;
    maxorder=20;
    z1=estimate_order_CGC( Y, kfold, maxorder,n_timepoints);
    pmax=z1.order; 
     
     
     %% estimation of CGC ----------------------------------------------1
    [CGC1,pCGC1]=CGC(Y,n_timepoints,pmax,[],0);
     CGC1= CGC1.*(ones(M,M)-eye(M,M));
    pCGC1=pCGC1.*(ones(M,M)-eye(M,M));
    
%% estimation of GLasso-CGC ----------------------------------------------3  
    lambda_GLasso_CGC=0.01;
    [~,GLasso_CGC1,pGLasso_CGC1]=GLasso_CGC(Y,n_timepoints,pmax,[],lambda_GLasso_CGC,0);
    GLasso_CGC1= GLasso_CGC1.*(ones(M,M)-eye(M,M));
pGLasso_CGC1=pGLasso_CGC1.*(ones(M,M)-eye(M,M));
    
%% %% estimation of mBTS-CGC -----------------------------------------2
[mBTSCGC1,pmBTSCGC1]=mBTSCGCImatrix(Y,n_timepoints,pmax,1);
mBTSCGC1=mBTSCGC1';
pmBTSCGC1=pmBTSCGC1';    
mBTSCGC1= mBTSCGC1.*(ones(M,M)-eye(M,M));
pmBTSCGC1=pmBTSCGC1.*(ones(M,M)-eye(M,M));    

  
%%  estimation of NCGC
 %---------------rbf center selection ----------cross validation-----------------   
   kfold=5; %must greater than one;
   maxCenter=20;
   [z, C]= estimate_numC_NCGC( Y, kfold, maxCenter,pmax ,n_timepoints);


    rbf_C=z.C; 
    num_C=size(rbf_C,1);   
   [NCGC1,pNCGC1,rbf_Center31,Phi_g31]=NCGC(Y,n_timepoints,pmax,num_C,rbf_C,1);%数据标准化了
      NCGC1= NCGC1.*(ones(M,M)-eye(M,M));
      pNCGC1=pNCGC1.*(ones(M,M)-eye(M,M));
%%  estimation of NENCGC

    rbf_C=z.C; 
    num_C=size(rbf_C,1);   
   [NENCGC1,pNENCGC1,rbf_Center32,ind3,Phi_g32,ecC1,ecC2]=NENCGC(Y,n_timepoints,pmax,num_C,rbf_C,1);%数据标准化了
         NENCGC1= NENCGC1.*(ones(M,M)-eye(M,M));
      pNENCGC1=pNENCGC1.*(ones(M,M)-eye(M,M));


%% COMPUTE AUROC & AUPR
%不设置pval截断，直接用推得的矩阵
% CGC1=CGC1.*(pCGC1<alpha);
% GLasso_CGC1=GLasso_CGC1.*(pGLasso_CGC1<alpha);
% mBTSCGC1=mBTSCGC1.*(pmBTSCGC1<alpha);
% NCGC1=NCGC1.*(pNCGC1<alpha);
% NENCGC1=NCGC1.*(pNENCGC1<alpha);


%CGC
holdfigures=0;
CGC1=(CGC1-min(min(CGC1)))/(max(max(CGC1))-min(min(CGC1)));
[CGC_AUROC1, ~ ,CGC_AUPR1] = AUCs(CGC1 , B);

%GLasso CGC
GLasso_CGC1=(GLasso_CGC1-min(min(GLasso_CGC1)))/(max(max(GLasso_CGC1))-min(min(GLasso_CGC1)));
[GLasso_CGC_AUROC1, ~ ,GLasso_CGC_AUPR1] = AUCs(GLasso_CGC1 , B); 

%mBTS-CGC
 mBTSCGC2=(mBTSCGC1-min(min(mBTSCGC1)))/(max(max(mBTSCGC1))-min(min(mBTSCGC1)));
[mBTSCGC_AUROC1, ~ ,mBTSCGC_AUPR1]= AUCs(mBTSCGC2 , B);


%NCGC
NCGC1=(NCGC1-min(min(NCGC1)))/(max(max(NCGC1))-min(min(NCGC1)));
[NCGC_AUROC1, ~ ,NCGC_AUPR1]= AUCs(NCGC1, B);

%NENCGC
NENCGC1=(NENCGC1-min(min(NENCGC1)))/(max(max(NENCGC1))-min(min(NENCGC1)));
[NENCGC_AUROC1, ~ ,NENCGC_AUPR1]= AUCs(NENCGC1, B);

CGC_AUROC(i)=CGC_AUROC1;
CGC_AUPR(i)=CGC_AUPR1;
GLasso_CGC_AUROC(i)= GLasso_CGC_AUROC1;
GLasso_CGC_AUPR(i)=GLasso_CGC_AUPR1;
mBTSCGC_AUROC(i)=mBTSCGC_AUROC1;
mBTSCGC_AUPR(i)=mBTSCGC_AUPR1;
NCGC_AUROC(i)=NCGC_AUROC1;
NCGC_AUPR(i)=NCGC_AUPR1;
NENCGC_AUROC(i)=NENCGC_AUROC1;
NENCGC_AUPR(i)=NENCGC_AUPR1;



CGCM(:,:,i)=CGC1; pCGCM(:,:,i)=pCGC1;  
GLasso_CGCM(:,:,i)=GLasso_CGC1; pGLasso_CGCM(:,:,i)=pGLasso_CGC1;
mBTSCGCM(:,:,i)=mBTSCGC1; pmBTSCGCM(:,:,i)=pmBTSCGC1;    
NCGCM(:,:,i)=NCGC1; pNCGCM(:,:,i)=pNCGC1;  
NENCGCM(:,:,i)=NENCGC1; pNENCGCM(:,:,i)=pNENCGC1;  


%% COMPUTE AUROC & AUPR
%不设置pval截断，直接用推得的矩阵
CGC2=CGC1.*(pCGC1<alpha);
GLasso_CGC2=GLasso_CGC1.*(pGLasso_CGC1<alpha);
mBTSCGC2=mBTSCGC1.*(pmBTSCGC1<alpha);
NCGC2=NCGC1.*(pNCGC1<alpha);
NENCGC2=NENCGC1.*(pNENCGC1<alpha);


[~,~,~,~,CGC_TPR(i), CGC_FPR(i), CGC_PPV(i), CGC_Fscore(i), CGC_ACC(i), CGC_MCC(i)] = cal_statistics(CGC2, B);
[~,~,~,~,GLasso_CGC_TPR(i),GLasso_CGC_FPR(i), GLasso_CGC_PPV(i), GLasso_CGC_Fscore(i), GLasso_CGC_ACC(i),GLasso_CGC_MCC(i)] = cal_statistics(GLasso_CGC2, B);
[~,~,~,~,mBTSCGC_TPR(i),mBTSCGC_FPR(i), mBTSCGC_PPV(i), mBTSCGC_Fscore(i), mBTSCGC_ACC(i), mBTSCGC_MCC(i)] = cal_statistics(mBTSCGC2, B);
[~,~,~,~,NCGC_TPR(i), NCGC_FPR(i), NCGC_PPV(i),NCGC_Fscore(i), NCGC_ACC(i), NCGC_MCC(i)] = cal_statistics(NCGC2, B);
[~,~,~,~,NENCGC_TPR(i), NENCGC_FPR(i), NENCGC_PPV(i), NENCGC_Fscore(i), NENCGC_ACC(i), NENCGC_MCC(i)] = cal_statistics(NENCGC2, B);


end

Aver_CGC_AUROC=sum(CGC_AUROC)/numsimu;   Std_CGC_AUROC =std(CGC_AUROC);
Aver_CGC_AUPR=sum(CGC_AUPR)/numsimu;     Std_CGC_AUPR =std(CGC_AUPR);
Aver_CGC_TPR=sum(CGC_TPR)/numsimu;       Std_CGC_TPR =std(CGC_TPR);
Aver_CGC_FPR=sum(CGC_FPR)/numsimu;       Std_CGC_FPR =std(CGC_FPR);
Aver_CGC_PPV=sum(CGC_PPV)/numsimu;       Std_CGC_PPV =std(CGC_PPV);
Aver_CGC_Fscore=sum(CGC_Fscore)/numsimu; Std_CGC_Fscore =std(CGC_Fscore);
Aver_CGC_ACC =sum(CGC_ACC)/numsimu;      Std_CGC_ACC =std(CGC_ACC);
Aver_CGC_MCC =sum(CGC_MCC)/numsimu;      Std_CGC_MCC=std(CGC_MCC);



Aver_GLasso_CGC_AUROC=sum( GLasso_CGC_AUROC)/numsimu;  Std_GLasso_CGC_AUROC=std(GLasso_CGC_AUROC);
Aver_GLasso_CGC_AUPR=sum(GLasso_CGC_AUPR)/numsimu;     Std_GLasso_CGC_AUPR=std(GLasso_CGC_AUPR);
Aver_GLasso_CGC_TPR=sum(GLasso_CGC_TPR)/numsimu;       Std_GLasso_CGC_TPR=std(GLasso_CGC_TPR);
Aver_GLasso_CGC_FPR=sum(GLasso_CGC_FPR)/numsimu;       Std_GLasso_CGC_FPR=std(GLasso_CGC_FPR);
Aver_GLasso_CGC_PPV=sum(GLasso_CGC_PPV)/numsimu;       Std_GLasso_CGC_PPV=std(GLasso_CGC_PPV);
Aver_GLasso_CGC_Fscore=sum(GLasso_CGC_Fscore)/numsimu; Std_GLasso_CGC_Fscore=std(GLasso_CGC_Fscore);
Aver_GLasso_CGC_ACC=sum(GLasso_CGC_ACC)/numsimu;       Std_GLasso_CGC_ACC=std(GLasso_CGC_ACC);
Aver_GLasso_CGC_MCC=sum(GLasso_CGC_MCC)/numsimu;       Std_GLasso_CGC_MCC=std(GLasso_CGC_MCC);

Aver_mBTSCGC_AUROC=sum(mBTSCGC_AUROC)/numsimu;    Std_mBTSCGC_AUROC=std(mBTSCGC_AUROC);
Aver_mBTSCGC_AUPR=sum(mBTSCGC_AUPR)/numsimu;      Std_mBTSCGC_AUPR=std(mBTSCGC_AUPR);
Aver_mBTSCGC_TPR=sum(mBTSCGC_TPR)/numsimu;        Std_mBTSCGC_TPR=std(mBTSCGC_TPR);
Aver_mBTSCGC_FPR=sum(mBTSCGC_FPR)/numsimu;        Std_mBTSCGC_FPR=std(mBTSCGC_FPR);
Aver_mBTSCGC_PPV=sum(mBTSCGC_PPV)/numsimu;        Std_mBTSCGC_PPV=std(mBTSCGC_PPV);
Aver_mBTSCGC_Fscore=sum(mBTSCGC_Fscore)/numsimu;  Std_mBTSCGC_Fscore=std(mBTSCGC_Fscore);
Aver_mBTSCGC_ACC=sum(mBTSCGC_ACC)/numsimu;        Std_mBTSCGC_ACC=std(mBTSCGC_ACC);
Aver_mBTSCGC_MCC=sum(mBTSCGC_MCC)/numsimu;        Std_mBTSCGC_MCC=std(mBTSCGC_MCC);

Aver_NCGC_AUROC=sum(NCGC_AUROC)/numsimu;         Std_NCGC_AUROC=std(NCGC_AUROC);
Aver_NCGC_AUPR=sum(NCGC_AUPR)/numsimu;           Std_NCGC_AUPR=std(NCGC_AUPR);
Aver_NCGC_TPR=sum(NCGC_TPR)/numsimu;             Std_NCGC_TPR=std(NCGC_TPR);
Aver_NCGC_FPR=sum(NCGC_FPR)/numsimu;             Std_NCGC_FPR=std(NCGC_FPR);
Aver_NCGC_PPV=sum(NCGC_PPV)/numsimu;             Std_NCGC_PPV=std(NCGC_PPV);
Aver_NCGC_Fscore=sum(NCGC_Fscore)/numsimu;       Std_NCGC_Fscore=std(NCGC_Fscore);
Aver_NCGC_ACC=sum(NCGC_ACC)/numsimu;             Std_NCGC_ACC=std(NCGC_ACC);
Aver_NCGC_MCC=sum(NCGC_MCC)/numsimu;             Std_NCGC_MCC=std(NCGC_MCC);

Aver_NENCGC_AUROC=sum(NENCGC_AUROC)/numsimu;      Std_NENCGC_AUROC=std(NENCGC_AUROC);
Aver_NENCGC_AUPR=sum(NENCGC_AUPR)/numsimu;        Std_NENCGC_AUPR=std(NENCGC_AUPR);
Aver_NENCGC_TPR=sum(NENCGC_TPR)/numsimu;          Std_NENCGC_TPR=std(NENCGC_TPR);
Aver_NENCGC_FPR=sum(NENCGC_FPR)/numsimu;          Std_NENCGC_FPR=std(NENCGC_FPR);
Aver_NENCGC_PPV=sum(NENCGC_PPV)/numsimu;          Std_NENCGC_PPV=std(NENCGC_PPV);
Aver_NENCGC_Fscore=sum(NENCGC_Fscore)/numsimu;    Std_NENCGC_Fscore=std(NENCGC_Fscore);
Aver_NENCGC_ACC=sum(NENCGC_ACC)/numsimu;          Std_NENCGC_ACC=std(NENCGC_ACC);
Aver_NENCGC_MCC=sum(NENCGC_MCC)/numsimu;          Std_NENCGC_MCC=std(NENCGC_MCC);


		
Result1=[Aver_CGC_AUROC            Std_CGC_AUROC;
         Aver_CGC_AUPR 	           Std_CGC_AUPR;
		 Aver_CGC_TPR              Std_CGC_TPR;
         Aver_CGC_FPR              Std_CGC_FPR;
         Aver_CGC_PPV              Std_CGC_PPV;
%         Aver_CGC_Fscore           Std_CGC_Fscore;
         Aver_CGC_ACC              Std_CGC_ACC;
%          Aver_CGC_MCC              Std_CGC_MCC;
         Aver_GLasso_CGC_AUROC     Std_GLasso_CGC_AUROC;
         Aver_GLasso_CGC_AUPR      Std_GLasso_CGC_AUPR;
         Aver_GLasso_CGC_TPR       Std_GLasso_CGC_TPR;
         Aver_GLasso_CGC_FPR       Std_GLasso_CGC_FPR;
         Aver_GLasso_CGC_PPV       Std_GLasso_CGC_PPV;
  %       Aver_GLasso_CGC_Fscore    Std_GLasso_CGC_Fscore;
         Aver_GLasso_CGC_ACC       Std_GLasso_CGC_ACC;
%          Aver_GLasso_CGC_MCC       Std_GLasso_CGC_MCC; 
         Aver_mBTSCGC_AUROC        Std_mBTSCGC_AUROC; 
         Aver_mBTSCGC_AUPR         Std_mBTSCGC_AUPR;
         Aver_mBTSCGC_TPR          Std_mBTSCGC_TPR;
         Aver_mBTSCGC_FPR          Std_mBTSCGC_FPR;
         Aver_mBTSCGC_PPV          Std_mBTSCGC_PPV;
%          Aver_mBTSCGC_Fscore       Std_mBTSCGC_Fscore; 
         Aver_mBTSCGC_ACC          Std_mBTSCGC_ACC;
%          Aver_mBTSCGC_MCC          Std_mBTSCGC_MCC;
         Aver_NCGC_AUROC           Std_NCGC_AUROC;
         Aver_NCGC_AUPR            Std_NCGC_AUPR;
         Aver_NCGC_TPR             Std_NCGC_TPR;
         Aver_NCGC_FPR             Std_NCGC_FPR;
         Aver_NCGC_PPV             Std_NCGC_PPV;
%          Aver_NCGC_Fscore          Std_NCGC_Fscore; 
         Aver_NCGC_ACC             Std_NCGC_ACC;
%          Aver_NCGC_MCC             Std_NCGC_MCC;
         Aver_NENCGC_AUROC         Std_NENCGC_AUROC; 
         Aver_NENCGC_AUPR          Std_NENCGC_AUPR; 
         Aver_NENCGC_TPR           Std_NENCGC_TPR; 
         Aver_NENCGC_FPR           Std_NENCGC_FPR; 
         Aver_NENCGC_PPV           Std_NENCGC_PPV; 
%          Aver_NENCGC_Fscore        Std_NENCGC_Fscore; 
         Aver_NENCGC_ACC           Std_NENCGC_ACC]; 
%          Aver_NENCGC_MCC           Std_NENCGC_MCC]; 
		 
Result=Result1';		 
		 
ind=0;
for i=1:M
    for j=1:M
        if i~=j
    ind=ind+1;
    pCGC_sum(i,j)=sum(pCGCM(i,j,:)<alpha);
    pGLasso_CGC_sum(i,j)=sum(pGLasso_CGCM(i,j,:)<alpha);
    pmBTSCGC_sum(i,j)=sum(pmBTSCGCM(i,j,:)<alpha);
    pNCGC_sum(i,j)=sum(pNCGCM(i,j,:)<alpha);
    pNENCGC_sum(i,j)=sum(pNENCGCM(i,j,:)<alpha);
        end
    end
end



%% plot true matrix
figure;
imagesc(B); 
set(gca,'XTickLabel',{}); 
set(gca,'YTickLabel',{});
title('Connectivity','FontName','Times New Roman','FontSize',14);

%% plot number of discovery rate matrix
figure;
subplot(151);
imagesc(pCGC_sum );
% set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10])% X坐标轴刻度数据点位置  
set(gca,'XTickLabel',{}) %X坐标轴刻度处显示的字符 
set(gca,'FontName','Times New Roman','FontSize',14)%设置坐标轴刻度字体名

% set(gca, 'YTick', [1 2 3 4 5 6 7 8 9 10])% X坐标轴刻度数据点位置  
set(gca,'YTickLabel',{}) %X坐标轴刻度处显示的字符 
set(gca,'FontName','Times New Roman','FontSize',14)%设置坐标轴刻度字体名
 title('CGC','FontName','Times New Roman','FontSize',14);

subplot(152);
imagesc(pGLasso_CGC_sum );
% set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10])% X坐标轴刻度数据点位置  
set(gca,'XTickLabel',{}) %X坐标轴刻度处显示的字符 
set(gca,'FontName','Times New Roman','FontSize',14)%设置坐标轴刻度字体名

% set(gca, 'YTick', [1 2 3 4 5 6 7 8 9 10])% X坐标轴刻度数据点位置  
set(gca,'YTickLabel',{}) %X坐标轴刻度处显示的字符 
set(gca,'FontName','Times New Roman','FontSize',14)%设置坐标轴刻度字体名
 title('GLasso-CGC','FontName','Times New Roman','FontSize',14);
 
subplot(153);
imagesc(pmBTSCGC_sum );
% set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10])% X坐标轴刻度数据点位置  
set(gca,'XTickLabel',{}) %X坐标轴刻度处显示的字符 
set(gca,'FontName','Times New Roman','FontSize',14)%设置坐标轴刻度字体名

% set(gca, 'YTick', [1 2 3 4 5 6 7 8 9 10])% X坐标轴刻度数据点位置  
set(gca,'YTickLabel',{}) %X坐标轴刻度处显示的字符 
set(gca,'FontName','Times New Roman','FontSize',14)%设置坐标轴刻度字体名
title('mBTS-CGC','FontName','Times New Roman','FontSize',14); 
 
 
subplot(154);
imagesc(pNCGC_sum);
% set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10])% X坐标轴刻度数据点位置  
set(gca,'XTickLabel',{}) %X坐标轴刻度处显示的字符 
set(gca,'FontName','Times New Roman','FontSize',14)%设置坐标轴刻度字体名

% set(gca, 'YTick', [1 2 3 4 5 6 7 8 9 10])% X坐标轴刻度数据点位置  
set(gca,'YTickLabel',{}) %X坐标轴刻度处显示的字符 
set(gca,'FontName','Times New Roman','FontSize',14)%设置坐标轴刻度字体名
title('NCGC','FontName','Times New Roman','FontSize',14);

subplot(155);
imagesc(pNENCGC_sum);
% set(gca, 'XTick', [1 2 3 4 5 6 7 8 9 10])% X坐标轴刻度数据点位置  
set(gca,'XTickLabel',{}) %X坐标轴刻度处显示的字符 
set(gca,'FontName','Times New Roman','FontSize',14)%设置坐标轴刻度字体名

% set(gca, 'YTick', [1 2 3 4 5 6 7 8 9 10])% X坐标轴刻度数据点位置  
set(gca,'YTickLabel',{}) %X坐标轴刻度处显示的字符 
set(gca,'FontName','Times New Roman','FontSize',14)%设置坐标轴刻度字体名
title('NENCGC','FontName','Times New Roman','FontSize',14);











