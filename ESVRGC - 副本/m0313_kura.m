%[TrainingTime, TestingTime, TrainingAccuracy, TestingAccuracy,Y,TY ]= esvr( [y2(1:3000,1),x2(1:3000,:)],[y2(3001:4090,1),x2(3001:4090,:)], 100, 'sig',10);

% load('C:\Users\Dell\Desktop\hww\SVR0927\SVR0927\SVR\x2.mat');
% load('C:\Users\Dell\Desktop\hww\SVR0927\SVR0927\SVR\y2.mat');

% load('C:\Users\Dell\Desktop\hww\SVR0927\SVR0927\SVR\simulation_difflen\matlab.mat');
% y = X1(1:50,:);

% load('C:\Users\Dell\Desktop\hww\SVR0927\SVR0927\SVR\simulation_difflen\10.6SMk3d0.01N1000_4.mat');
% y = XX';

% load('E:\hww\SVR0927\SVR0927\SVR\simulation_difflen\realization_linear_4096.mat');
% data = XX;



% load('E:\hww\ESVRGC\simulation_difflen\10.6SMk3d0.03N1000_5.mat'); 
load('E:\hww\ESVRGC\simulation_difflen\0315.mat'); 
y = XX';

% load('E:\hww\ESVRGC\simulation_difflen\matlab.mat');
% y = data';
pmax=1;
D=pmax+1;
n_time_points = 2;
[y_S,~,Phi_g3,ind3,C,ecC1,ecC2]=data_process_mixedembedding(y,n_time_points,D,1,[],1);
N=size(y,2);
m=size(y,1)-pmax;
tmp=zeros(m,pmax,N);
tmp1=cell(N,1);
            

   


 stdData = zscore(y);
MinMaxData = mapminmax(stdData');
MinMaxData = MinMaxData';

sequence_length = pmax;
num_shift = 2;
x2 = zeros(size(y,1)-sequence_length,sequence_length,size(y,2));
y2 = zeros(size(y,1)-sequence_length,size(y,2));
for p =1:num_shift:size(MinMaxData,1)
    
    if (p+ sequence_length  >size(MinMaxData,1))
        break;
    end
    x2(p,:,:) = [MinMaxData(p:p + sequence_length-1,:)];
    
    y2(p,:) = [MinMaxData(p + sequence_length,:)];
%     x2 = [x2, y(p:p + sequence_length-1,:)];
%     
%     y2 = [y2, y(p + sequence_length,:)];
end
    
    
    
    
    
    for j =1:N
        tmp=zeros(m,sequence_length,N);
    for i =1:m
        for k =1:N
        idx = find(ecC2{j}(:,1)==k);
        if~isempty(idx)
            klags = ecC2{j}(idx,2);
            tmp(i,klags,k)= x2(i,klags,k);
        end
        end
    end
    tmp1{j}=tmp;
    end

    
    hiddenneus = 50;
tolerC = 10000;

granger_matrix = zeros(size(x2,3),size(x2,3));
var_denominator = zeros(1,size(x2,3));
trains = round(0.9*size(x2,1));
tests = size(x2,1)-trains;
for j = 1:1:size(y,2)
    x2 = tmp1{j};
    tmp_y =  y2(1:trains,j);
    tmp_y_tests = y2(trains+1:size(x2,1),j);
    tmp_x =  x2(1:trains,:,:);
    tmp_x_tests = x2(trains+1:size(x2,1),:,:);
    input_set = [];
    for n = 1:1:size(y,2)
        input_set = [input_set,n];
    end
    x_11 = tmp_x(:,:,input_set);
    x_111 = reshape(x_11,[size(tmp_x,1),sequence_length * size(y,2)]);
    x_11_test = tmp_x_tests(:,:,input_set);
    x_111_test = reshape(x_11_test,[size(tmp_x_tests,1),sequence_length*size(y,2)]);
    [TrainingTime, TestingTime, TrainingAccuracy, TestingAccuracy,TY] = esvrkernel( [tmp_y,x_111],[tmp_y_tests,x_111_test],tolerC, 'RBF_kernel',10);
    
    var_denominator(1,j) = TrainingAccuracy;
    
    for k = 1:1:size(y,2)
        x_11 = tmp_x(:,:,input_set);
        channel_del_idx = input_set(k);
        x_11(:,:,channel_del_idx) = 0;
        x_111 = reshape(x_11,[size(tmp_x,1),sequence_length*size(y,2)]);
        x_11_test = tmp_x_tests(:,:,input_set);
        x_11_test(:,:,channel_del_idx) = 0;
        x_111_test = reshape(x_11_test,[size(tmp_x_tests,1),sequence_length * size(y,2)]);
        [TrainingTime_1, TestingTime_1, TrainingAccuracy_1, TestingAccuracy_1,TY_1] = esvrkernel( [tmp_y,x_111],[tmp_y_tests,x_111_test], tolerC,'RBF_kernel', 10);
        granger_matrix(k,j) = TrainingAccuracy_1;
        
    end
end
    
granger_matrix = granger_matrix ./ repmat(var_denominator,size(var_denominator,2),1);
granger_matrix(granger_matrix < 1 ) =1;
for i =1:1:size(y,2)
    granger_matrix(i,i) = 1;
end

granger_matrix = log(granger_matrix);
% granger_matrix=(granger_matrix-min(min(granger_matrix)))/(max(max(granger_matrix))-min(min(granger_matrix)));

% HeatMap(flipud(granger_matrix));
imagesc(granger_matrix');
set(gcf,'colormap',jet);
caxis([0 max(max(granger_matrix))]);
colorbar;
% figure;
