function xM = SysSchelter06bm1(n)
% xM = SysSchelter06bm1(n)
% SysSchelter06bm1 generates a multivariate linear time series of 5 
% variables from a VAR(4) model. In parenthesis the largest lag of the
% driving variable in the expression for the response variable.
% x1(4)-->x2  x1(2)-->x4  
% x2(2)-->x4 
% x4(1)-->x5
% x5(1)-->x1  x5(2)-->x2  x5(3)-->x3
% INPUTS
% - n   : the time series length
% OUTPUTS
% - xM  : the n x 5 matrix of the generated time series
% Reference 
% Schelter, Winterhalder, Hellwig, Guschlbauer, L¨¹cking, Timmer
% Direct or indirect? Graphical models for neural oscillators
% J Physiology - Paris 99:37-46, 2006.
% Example Model 1: Five-dimensional vector autoregressive process of order
% 4
ntrans = 100; % transient period
K = 5; % number of variables
P = 4; % order of VAR
wM = randn(n+ntrans,K);
xM = NaN*ones(n+ntrans,K);
xM(1:P,:) = wM(1:P,:);
for t=P+1:n+ntrans
    xM(t,1)=0.4*xM(t-1,1)-0.5*xM(t-2,1)+0.4*xM(t-1,5)+wM(t,1);
    xM(t,2)=0.4*xM(t-1,2)-0.3*xM(t-4,1)+0.4*xM(t-2,5)+wM(t,2);
    xM(t,3)=0.5*xM(t-1,3)-0.7*xM(t-2,3)-0.3*xM(t-3,5)+wM(t,3);
    xM(t,4)=0.8*xM(t-3,4)+0.4*xM(t-2,1)+0.3*xM(t-2,2)+wM(t,4);
    xM(t,5)=0.7*xM(t-1,5)-0.5*xM(t-2,5)-0.4*xM(t-1,4)+wM(t,5);
end;
xM = xM(ntrans+1:n+ntrans,:);