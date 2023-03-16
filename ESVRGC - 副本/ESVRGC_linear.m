%[TrainingTime, TestingTime, TrainingAccuracy, TestingAccuracy,Y,TY ]= esvr( [y2(1:3000,1),x2(1:3000,:)],[y2(3001:4090,1),x2(3001:4090,:)], 100, 'sig',10);

% load('C:\Users\Dell\Desktop\hww\SVR0927\SVR0927\SVR\x2.mat');
% load('C:\Users\Dell\Desktop\hww\SVR0927\SVR0927\SVR\y2.mat');

% load('C:\Users\Dell\Desktop\hww\SVR0927\SVR0927\SVR\simulation_difflen\matlab.mat');
% y = X1(1:50,:);

% load('C:\Users\Dell\Desktop\hww\SVR0927\SVR0927\SVR\simulation_difflen\10.6SMk3d0.01N1000_4.mat');
% y = XX';

load('E:\hww\SVR0927\SVR0927\SVR\simulation_difflen\realization_linear_4096.mat');
y = data';

% 
[y_S,~,Phi_g3,ind3,C,ecC1,ecC2]=data_process_mixedembedding(y,4096,4,6,[],1);

% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% 
% %load('C:\Users\Dell\Desktop\hww\SVR0927\SVR0927\SVR\y2.mat');
% stdData = zscore(y);
% MinMaxData = mapminmax(stdData');
% MinMaxData = MinMaxData';
% 
% sequence_length = 5;
% num_shift = 1;

for k = 1:1:size(y,2)
    y1 = y(:,ind3{k});
stdData = zscore(y1);
MinMaxData = mapminmax(stdData');
MinMaxData = MinMaxData';

sequence_length = 5;
num_shift = 1;
    x2 = zeros(size(y,1)-sequence_length,sequence_length,ind3{k});
    y2 = zeros(size(y,1)-sequence_length,ind3{k});
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
    %x2 = reshape(x2,[size(y,1)-sequence_length,sequence_length,size(y,2)]);
    %y2 = reshape(y2,[size(y,1)-sequence_length,sequence_length]);

    hiddenneus = 50;
    tolerC = 1;

    granger_matrix = zeros(size(x2,3),size(x2,3));
    var_denominator = zeros(1,size(x2,3));
    trains = round(0.9*size(x2,1));
    tests = size(x2,1)-trains;
      tmp_y =  y2(1:trains,ind3{k});
      tmp_y_tests = y2(trains+1:size(x2,1),ind3{k});
      tmp_x =  x2(1:trains,:,ind3{k});
      tmp_x_tests = x2(trains+1:size(x2,1),:,ind3{k});
%       input_set = [];
%       for i = 1:1:size(x2,3)
%         input_set = [input_set,i];
%       end
      x_11 = tmp_x(:,:,[ind3{k}]);
      x_111 = reshape(x_11,[size(tmp_x,1),size(ind3{k},2)*size(y,2)]);
      x_11_test = tmp_x_tests(:,:,ind3{k});
      x_111_test = reshape(x_11_test,[size(tmp_x_tests,1),size(ind3{k},2)*size(y,2)]);
      
      
      
%     tmp_y =  y2(1:trains,k);
%     tmp_y_tests = y2(trains+1:size(x2,1),k);
%     tmp_x =  x2(1:trains,:,:);
%     tmp_x_tests = x2(trains+1:size(x2,1),:,:);
%     input_set = [];
%     for i = 1:1:size(x2,3)
%         input_set = [input_set,i];
%     end
%     x_11 = tmp_x(:,:,input_set);
%     x_111 = reshape(x_11,[size(tmp_x,1),sequence_length*size(y,2)]);
%     x_11_test = tmp_x_tests(:,:,input_set);
%     x_111_test = reshape(x_11_test,[size(tmp_x_tests,1),sequence_length*size(y,2)]);
    [TrainingTime, TestingTime, TrainingAccuracy, TestingAccuracy,Y,TY] = esvr( [tmp_y,x_111],[tmp_y_tests,x_111_test], hiddenneus, 'sig',tolerC);
    
    var_denominator(1,k) = TrainingAccuracy;
    
    for j = 1:1:ind3{k}
        x_11 = tmp_x(:,:,ind3{k});
        if (ind3{k}(j)==1)
            break;
        end
           
        channel_del_idx = ind3{k}(j);
        x_11(:,:,channel_del_idx) = 0;
        x_111 = reshape(x_11,[size(tmp_x,1),size(ind3{k},2)*size(y1,2)]);
        x_11_test = tmp_x_tests(:,:,ind3{k});
        x_11_test(:,:,channel_del_idx) = 0;
        x_111_test = reshape(x_11_test,[size(tmp_x_tests,1),size(ind3{k},2)*size(y,2)]);
        [TrainingTime_1, TestingTime_1, TrainingAccuracy_1, TestingAccuracy_1,Y_1,TY_1] = esvr( [tmp_y,x_111],[tmp_y_tests,x_111_test], hiddenneus, 'sig',tolerC);
        granger_matrix(j,k) = TrainingAccuracy_1;
    end
end
granger_matrix = granger_matrix ./ repmat(var_denominator,size(var_denominator,2),1);
granger_matrix(granger_matrix < 1 ) =1;
for i =1:1:size(y,2)
    granger_matrix(i,i) = 1;
end

granger_matrix = log(granger_matrix);
HeatMap(flipud(granger_matrix));
%HeatMap(flipud(A));