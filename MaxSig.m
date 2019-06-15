%% Find Max Signature

clear all; close all; clc;

%% Open file 

in=input('Name of the file: ');
filename=strcat(in,'_sig.xlsx');

[a,b,raw] = xlsread(filename);


%% Split data

signature=raw(2:end,1);
patient=raw(1,2:end);
values=raw(2:end,2:end);

%% Processing

pat_sig1=[];
pat_sig2=[];
pat_sig3=[];
pat_sig4=[];
pat_sig5=[];

for i=1:length(patient)
    [val, idx] = max(cell2mat(values(:,i)));
    if idx==1
        pat_sig1=[pat_sig1 patient(i)];
    elseif idx==2
        pat_sig2=[pat_sig2 patient(i)];
    elseif idx==3
        pat_sig3=[pat_sig3 patient(i)];
    elseif idx==4
        pat_sig4=[pat_sig4 patient(i)];
    elseif idx==5
        pat_sig5=[pat_sig5 patient(i)];
    end
end

%% Prepare Table Output
lp1=length(pat_sig1);
lp2=length(pat_sig2);
lp3=length(pat_sig3);
lp4=length(pat_sig4);
lp5=length(pat_sig5);
maxl=max([lp1,lp2,lp3,lp4,lp5]);
output=cell(maxl+1,length(signature));
output(1,:)=signature;
output(2:lp1+1,1)=pat_sig1;
output(2:lp2+1,2)=pat_sig2;
output(2:lp3+1,3)=pat_sig3;
output(2:lp4+1,4)=pat_sig4;
output(2:lp5+1,5)=pat_sig5;
%% Send to excel

%outputname=strcat(in,'MaxSig.xlsx');

%xlswrite(outputname,output);

%% Make Bar Plots

c = categorical(signature);
l_s=length(signature);
if l_s==3
    bar(c,[lp1,lp2,lp3])
elseif l_s==4
    bar(c,[lp1,lp2,lp3,lp4])
elseif l_s==5
    bar(c,[lp1,lp2,lp3,lp4,lp5])
end
ylabel('# of patients with signature')
title(strcat(in,' patients providing to the signature'))
