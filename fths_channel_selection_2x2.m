% This program realizes function of channel selection (two-by-two pairs).

clear;
clc;

% Define output matrix
output{1,1}={'Positive Channel 1'};
output{1,2}={'Positive Channel 2'};
output{1,3}={'Negative Channel 1'};
output{1,4}={'Negative Channel 2'};
output{1,5}={'Accuracy'};

% Define counter of output matrix
K=2;

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
    
    % Generate all channel combinations
    combination_num=2;
    channel_c=combntns(1:62,combination_num*2);
    [row_c,~]=size(channel_c);
    
    
    
    
    % Calculate the classification results under various channel combination conditions in turn
    for m=1:row_c
        
        % Perform a second combination of the current combination number
        % to determine all possible combinations of the positive channels
        channel_pos=combntns(channel_c(m,:),combination_num);
        [row_p,~]=size(channel_pos);
        
        for n=1:row_p
            
            % Calculate the corresponding negative channels according to the
            % current positive channels
            channel_neg=channel_c(m,:);
            [~,ma,~]=intersect(channel_neg,channel_pos(n,:));
            channel_neg(ma)=[];
            
            for time=1:column
                hs(1,1)=abs(eegdata(channel_pos(n,1),time)-positive_center(1));
                hs(1,2)=abs(eegdata(channel_pos(n,1),time)-positive_center(2));
                ......
                    hs(1,i)=abs(eegdata(channel_pos(n,1),time)-positive_center(i));
                
                hs(2,1)=abs(eegdata(channel_pos(n,2),time)-positive_center(1));
                hs(2,2)=abs(eegdata(channel_pos(n,2),time)-positive_center(2));
                ......
                    hs(2,i)=abs(eegdata(channel_pos(n,2),time)-positive_center(i));
                
                ss(1,1)=abs(eegdata(channel_neg(1),time)-negative_center(1));
                ss(1,2)=abs(eegdata(channel_neg(1),time)-negative_center(2));
                ......
                    ss(1,j)=abs(eegdata(channel_neg(1),time)-negative_center(j));
                
                ss(2,1)=abs(eegdata(channel_neg(2),time)-negative_center(1));
                ss(2,2)=abs(eegdata(channel_neg(2),time)-negative_center(2));
                ......
                    ss(2,j)=abs(eegdata(channel_neg(2),time)-negative_center(j));
                
                % Using Euclidean distance to judge the emotional state of current segment
                if sum(sum(hs))<sum(sum(ss))
                    current=1;
                else current=0;
                end;
                
                % Compare the current results with the state in which the eeg data should be located.
                % If they are the same, the correct counter +1, otherwise the error counter +1
                if ~xor(current,e)
                    result(time)=1;
                else result(time)=0;
                end;
                
            end;
            
            % Record result
            output{K,1}=channel_pos(n,1);
            output{K,2}=channel_pos(n,2);
            output{K,3}=channel_neg(1);
            output{K,4}=channel_neg(2);
            output{K,5}=sum(result)/time*100;
            K=K+1;
            
        end;
    end;
end;
