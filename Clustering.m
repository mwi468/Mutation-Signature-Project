%% Clustering

clear all; close all; clc;

%% Prepare data

in=input('Name of the file: ');
filename=strcat(in,'TNM.xlsx');

[a,b,raw] = xlsread(filename);

%%
groups=raw(1,2:97)';
patient=raw(2:end,1);
matrix=cell2mat(raw(2:end,2:97));

%% Distance
distance=pdist(matrix);

link=linkage(distance);
T=cluster(link,'maxclust',4);
cutoff=median([link(end-2,3) link(end-1,3)]);
figure()
dendrogram(link,4,'ColorThreshold',cutoff)
title(strcat(in, ' patients clustering'))

%% Organize Cluster of IDs

r1={};
r2={};
r3={};
r4={};
for i=1:length(T)
    if T(i)==1
        r1=[r1,patient(i)];
    elseif T(i)==2
        r2=[r2,patient(i)];
    elseif T(i)==3
        r3=[r3,patient(i)];
    elseif T(i)==4
        r4=[r4,patient(i)];
    end
end

out = cell(length(T),4);
for row = 1 : length(r1)
  out{row,1} = r1{row};
end
for row = 1 : length(r2)
  out{row,2} = r2{row};
end
for row = 1 : length(r3)
  out{row,3} = r3{row};
end
for row = 1 : length(r4)
  out{row,4} = r4{row};
end

%% Export to Excel
xlswrite(strcat(in,'_cluster.xlsx'),out)


