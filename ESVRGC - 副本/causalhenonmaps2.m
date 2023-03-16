function xM = causalhenonmaps2(m,c,n,aV,x0M)
% xM = causalhenonmaps2(c,n,aV,x0V)
% CAUSALHENONMAPS2 generates m time series from identical uni-directional
% coupled Henon maps with coupling strength c. Each map causes the next and 
% the last map does not cause any other map.
% INPUTS
% - m   : the number of coupled maps (time series)
% - c   : the coupling strength
% - n   : the time series length
% - aV  : the two parameters in a vector (if omitted the standard
%         parameters for chaos are chosen)
% - x0M : initial conditions, a matrix of size 2 x m
% OUTPUTS
% - xM  : the n x m matrix of the generated time series

if nargin==4
    x0M = NaN;
elseif nargin==3
    x0M = NaN;
    aV = [1.4 0.3];
end
if isempty(aV)
    aV = [1.4 0.3];
end  
if length(c)==1
    cV = c*ones(m,1);
else
    cV = c;
end
if isnan(x0M)
    ntrans = 1000;
    x0limV = [-1 1 -0.3 0.3]';
    xM = NaN*ones(n+ntrans,m);
    tmpM = rand(2,m);
    xM(1,:)=(x0limV(4)-tmpM(1,:))/(x0limV(4)-x0limV(3));
    xM(2,:)=(x0limV(2)-tmpM(2,:))/(x0limV(2)-x0limV(1));
else
    ntrans = 2;
    xM(1:2,:)= x0M;
end
for t=3:n+ntrans
    xM(t,1)= aV(1) - xM(t-1,1)^2 + aV(2)*xM(t-2,1);
    for im=2:m
        xM(t,im)= aV(1) - cV(im)*xM(t-1,im-1)*xM(t-1,im) - ...
            (1-cV(im))*xM(t-1,im)^2 + aV(2)*xM(t-2,im);
    end
end
xM = xM(ntrans+1:n+ntrans,:);
