% This program classifies the 1/f feature angles of all eeg data.
% knn and svm classifiers are used respectively.
clear;
clc;

% Define classification parameters
% G1 corresponds to knn, and G2 corresponds to svm
% Assignment is enabled at 1 and disabled at 0.
G1=0;
G2=0;

% ************************ knn classifcation *************************
if G1
    
    % Replace xxxx with your .mat or .dat file
    temp_data=load('xxxx.mat');
    temp_name=fieldnames(temp_data);
    
    for n=1:length(tempname)
        
        eegdata=temp.(tempname{n});
        
        if strfind(tempname{n},'positive')
            e=1;
        elseif strfind(tempname{n},'negative')
            e=0;
        else e=-1;
        end;
        
        
        % Suppose that the matrix storing EEG data is called eegdata
        [row,column]=size(eegdata);
        
        % Define knn clustering centers of positive and negative stimuli,
        % where pc? and nc? are the average feature angles of their stimuli.
        positive_center=[pc1 pc2 ... pci];
        negative_center=[nc1 nc2 ... ncj];
        
        result{1,1,n}={'Number of segments'};
        result{1,2,n}={'Number of correct segments'};
        result{1,3,n}={'Accuracy'};
        
        % 计算各组eeg数据的分类结果
        for j=1:column
            for i=1:row
                
                % 计算当前分段、当前导联1/f特征角度与正向/负向聚类中心的欧氏距离
                hr(i,j)=sum(abs(eegdata(i,j)-positive_center));
                sr(i,j)=sum(abs(eegdata(i,j)-negative_center));
                
                % 将当前分段、当前导联判断为距离小的状态
                if hr(i,j)<sr(i,j)
                    seg_para(i,j)=1;
                else seg_para(i,j)=0;
                end;
                
            end;
            
            % Calculate the statistical results of the current segment in all electrodes
            column_para=sum(seg_para(:,j))/row;
            
            % The current lead classification result is given according to the emotion parameter.
            switch e
                case 1
                    if column_para>=0.5
                        result_c(j)=1;
                    else result_c(j)=0;
                    end;
                case 0
                    if column_para<=0.5
                        result_c(j)=1;
                    else result_c(j)=0;
                    end;
            end;
        end;
        
        
        % Record results
        result{2,1}=column;
        result{2,2}=sum(result_c);
        result{2,3}=sum(result_c)/column*100;
    end;
end;
% ****************************************************************

% ------------------------ svm classification ----------------------------
if G2
    
    % Replace xxxx with your .mat or .dat file
    load xxxx.mat; % positive data
    load yyyy.mat; % negative data
    
    % Create positive/negative data matrices and emotional tags（positive is 1 and negative is 2）
    data_p=eegdata_positive; % Replace eegdata_positive with your eegdata
    [~,dim_p]=size(data_p);
    label_p=ones(dim_p,1);
    
    data_n=eegdata_negative; % Replace eegdata_negative with your eegdata
    [~,dim_n]=size(data_n);
    label_n=2*ones(dim_n,1);
    
    % Define the matrix of classification results
    output=zeros(2,2);
    
    % Define the final result matrix
    resultmatrix{1,1}='SVM';
    resultmatrix{1,2}='positive_out';
    resultmatrix{1,3}='negative_out';
    resultmatrix{2,1}='positive_in';
    resultmatrix{3,1}='negative_in';
    
    % Leave-one-out method
    for N=1:(dim_p+dim_n)
        
        tempdata1(:,:)=data_p(:,:);
        tempdata2(:,:)=data_n(:,:);
        train_data=[tempdata1';tempdata2'];
        test_data=train_data(N,:);
        train_data(N,:)=[];
        
        train_label=[label_p;label_n];
        test_label=train_label(N);
        train_label(N)=[];
        
        Factor = svmtrain(train_data, train_label);
        
        predict_label(N) = svmclassify(Factor, test_data);
        
        if test_label==1&&predict_label(N)==1
            output(1,1)=output(1,1)+1;
        elseif test_label==1&&predict_label(N)==2
            output(1,2)=output(1,2)+1;
        elseif test_label==2&&predict_label(N)==1
            output(2,1)=output(2,1)+1;
        else output(2,2)=output(2,2)+1;
        end;
        
    end;
    
    resultmatrix{2,2}=output(1,1)/sum(output(1,:))*100;
    resultmatrix{2,3}=output(1,2)/sum(output(1,:))*100;
    resultmatrix{3,2}=output(2,1)/sum(output(2,:))*100;
    resultmatrix{3,3}=output(2,2)/sum(output(2,:))*100;
end;
% ---------------------------------------------------------------------