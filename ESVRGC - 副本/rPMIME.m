% Example: Compute PMIME and PMIMEsig on a multiple time series generated
% by the coupled Henon maps
clear all
K = 5; % When data are to be generated, number of variables
C= 0.2;
n = 500; % The time series segment length
nnei = 5; % Number of nearest neighbors for PMIME
Lmax = 5; % The window size to search for components for MIME
thresh = 0.03; % The fixed hresholds for terminartion of MIME,PMIME (as proportion)
T = 1; % The prediction time for PMIME
nsur = 100;
alpha = 0.05; % The adjusted threshold for for termination of MIME,PMIME
toplot = 1;
cax = [0 0.07];
figno=0;

xM = coupledhenonmaps2(K,C,n);
maxite = 10;
itenow = 1;
while max(max(abs(xM))) > 5 && itenow <=maxite
    xM = coupledhenonmaps2(K,C,n);
end
if itenow>maxite
    error(sprintf('Could not generate valid coupled Henon data after %d trials.',maxite));
end

tstart1 = tic;
[RM,ecC] = PMIME(xM,Lmax,T,nnei,thresh,1);
telapsed1 = toc(tstart1);
for iK=1:K
    fprintf('===== i=%d \n',iK);
    disp(ecC{iK});
end
tstart2 = tic;
[R2M,ec2C] = PMIMEsig(xM,Lmax,T,nnei,nsur,alpha,1);
telapsed2 = toc(tstart2);
for iK=1:K
    fprintf('===== i=%d \n',iK);
    disp(ec2C{iK});
end
fprintf('Time PMIMEsig=%2.5f  Time PMIME=%2.5f \n',telapsed1,telapsed2);
fprintf('R computed by PMIME: \n')
disp(RM)
fprintf('R computed by PMIMEsig: \n')
disp(R2M)
