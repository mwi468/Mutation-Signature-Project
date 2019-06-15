%% 
clear all; close all; clc;
%% Open Excel Files

in=input('Name of the file: ');
filename=strcat(in,'.xlsx');

[a,b,raw] = xlsread(filename);

%% Process Data

for i=1:length(raw)
    if i==1
    elseif i==2
        patient=raw(2,1);
        pat_cnt=1;
        pat_mut=[1];
    elseif sum(strcmp(raw{i,1},patient))
        pat_mut(strcmp(raw{i,1},patient))=pat_mut(strcmp(raw{i,1},patient))+1;
    else
        pat_cnt=pat_cnt+1;
        pat_mut= [pat_mut 1];
        patient=[patient raw(i,1)];
    end
end

%% Display
data= [patient', num2cell(pat_mut)'];
Min= min(pat_mut)
Max= max(pat_mut)
Median= (Max+Min)/2
Avg=mean(pat_mut)
Std=std(pat_mut)
pat_cnt
total_mut=sum(pat_mut)

%% Organize Greater
greater=[];
for i=1:length(data)
    if cell2mat(data(i,2))>Median
        greater=[greater; data(i,:)];
    else
    end
end


%% Save Data

outputname=strcat(in,'_mut.xlsx');

xlswrite(outputname,greater(:,1));


        