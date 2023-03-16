% 
%   lambda_vec=[0.001,0.01, 0.1];
%    n_time_points = 20; 
%    D=2;%time-lags, D-1 is the model order
%    nonlinear=true;
%    maxCenter=5;
   
function [z ,C]= estimate_numC_NCGC( data, kfold, maxCenter ,pmax,n_time_points)
%==========================================================================
    N=size(data,2);
    D=pmax+1;%time-lags, D-1 is the model order
  
%     [~,N_t] = size(Phi_g2);
      

MSE=zeros(maxCenter,1);   
    C=cell(maxCenter,1);
for num_C= 1:maxCenter


    nonlinear=true;
    rbf_C=[];

    [y_S,~,Phi_g3,C1]=data_process(data,n_time_points,D,num_C,rbf_C,nonlinear);
     C{num_C}=C1;
    MSEerror = zeros(N,1);
         for i=1:N
             Y=y_S(:,i);
             

%% 
%              [~,N_t] = size(Phi_g2);
%              pool_Si =((i-1)*(D-1)+1):(i*(D-1));
%              ind_pool = sort(setdiff(1:N_t,pool_Si));      % The pool excluding columns corresponding to gene S(i)
%              Phi = Phi_g2(:,ind_pool);


             Phi = Phi_g3;
             [n, p] = size(Phi);
             CVO = cvpartition( n, 'kfold' , kfold );            
             err = zeros(CVO.NumTestSets,1);
 %%      
           for j = 1:CVO.NumTestSets
          
            % get ith training set
            trIdx = CVO.training(j);
          
          % get ith test set
            teIdx = CVO.test(j);
          
          % train Regression using training set
% 		      val = lasso(Phi(trIdx,:), Y(trIdx), 'Alpha', 1, 'Lambda',lambda_vec(m));
%                cardz = numel(find(val));

           [~, ~,Am]= Regression_LS_1( Y(trIdx),Phi(trIdx,:));
           val=Am';
           


          % test using testing set          
            Ypred = Phi(teIdx,:) * val;
                    
          % calculate MSE error            
            temp = Y(teIdx) - Ypred;
                    
          % clear ypred
            clear Ypred;
          
          % calculate average error for the ith test set
            err(j) = (temp'*temp)/length(temp);  
            %%%%%%%%%%%%%%
%              err(j) = (temp'*temp)/((1 - (cardz / (N-1)))^2) / (N-1); 
%                err(j) = (temp'*temp)/((1 - (cardz / N))^2) / N; 
          end
          % calculate the mean error over all test sets
            MSEerror(i) = mean(err);
         end  
       MSE(num_C)=mean(MSEerror);
        
end
    k=find(MSE==min(MSE));
    z.num_C=k; 
    z.mse=MSE;
%     k
%     C{k}
%     MSE
    z.C=C{k};%k每次是一个值 不能是2 3 4 5这么多个 把下面的end放到前面去  这个是做完循环之后才寻找最小的一个的
%end%
end


 function [Sigma, Up,Am]= Regression_LS_1(Y,Phi) % Linear Regression
 %% 
 %sigma : 误差方差
 %Up    ：回归误差
 %%
  if isempty(Phi) %if no conditioning, ce will be the Entropy of B
    Sigma=var(Y);
    Up=Y;
 else
    Yb=Y'; % inversion works with data organized in rows
    Z=Phi';
    Am=Yb/Z; % least squares!

    Yp=Am*Z; 
    Up=Yb-Yp;
    Sigma=cov(Up');
  end

 end

%
%==========================================================================
%               
%               PURPOSE:
%
%   Estimate L1 regularization parameter in LASSO via cross validation.
%
%   Algorithm for solving the Lasso problem:
%         0.5 * (y - X*beta)'*(y - X*beta) + lambda * ||beta||_1
%                                              
%   where ||beta||_1 is the L_1 norm i.e., ||beta||_1 = sum(abs( beta ))
%
%   We use the method proposed by Fu et. al based on single co-ordinate
%   descent. For more details see GP's notes or the following paper:
%
%   Penalized Regressions: The Bridge Versus the Lasso
%   Wenjiang J. FU, Journal of Computational and Graphical Statistics, 
%   Volume 7, Number 3, Pages 397?416, 1998
%   
%==========================================================================
%               
%               USAGE:
%
%               z = estimateLassoLambda( y, X, kfold, lambda_vec )
%
%==========================================================================
%
%==========================================================================
%               
%               INPUTS:
%
%       =>          y = n by 1 response vector
%
%       =>          X = n by p design matrix
%
%       =>      kfold = use kfold cross validation (e.g., kfold = 5)
%
%       => lambda_vec = a column vector of potential regularization
%                       parameters to consider
%
%==========================================================================
%
%==========================================================================
%               
%               OUTPUTS:
%       =>     z.lambda = optimal lambda value computed by kfold cross-validation
%    
%       => z.lambda_vec = user supplied vector of potential lambda values
%        
%       =>   z.MSEerror = average MSE over kfold folds for each lambda value
%  
%       =>  z.min_index = index of minimum of MSEerror
% 
%       =>        z.CVO = kfold cross-validation structure
%   
%       =>          z.y = user supplied response vectors
%    
%       =>          z.X = user supplied design matrix
% 
%       =>      z.kfold = user supplied value for kfold 
%
%==========================================================================
%
%==========================================================================
%   
%       Copyright 2011 : Gautam V. Pendse
%
%               E-mail : gautam.pendse@gmail.com
%
%                  URL : http://www.gautampendse.com
%
%==========================================================================