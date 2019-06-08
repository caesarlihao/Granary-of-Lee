% This program realizes function of channel selection (one-by-one pairs).

clear;
clc;

% Define output matrix
output{1,1}={'Positive Channel'};
output{1,2}={'Negative Channel'};
output{1,3}={'Accuracy'};

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
    
    % Select each channel as a positive emotional feature channel in turn,
    % and match one of the most suitable negative emotional feature channels from the remaining channels.
    for positive=1:row
        for negative=1:row
            
            if negative~=positive
                
                for time=1:column
                    hs(1)=abs(eegdata(positive,time)-positive_center(1));
                    hs(2)=abs(eegdata(positive,time)-positive_center(2));
                    ......
                        hs(i)=abs(eegdata(positive,time)-positive_center(i));
                    
                    ss(1)=abs(eegdata(positive,time)-negative_center(1));
                    ss(2)=abs(eegdata(positive,time)-negative_center(2));
                    ......
                        ss(j)=abs(eegdata(positive,time)-negative_center(j));
                    
                    % Using Euclidean distance to judge the emotional state of current segment
                    if sum(hs)<sum(ss)
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
                output{K,1}=positive;
                output{K,2}=negative;
                output{K,3}=sum(result)/time*100;
                K=K+1;
                
            end;
        end;
    end;
end;