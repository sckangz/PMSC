     f=load('bbc_seg14of4.mat');%=load('C:\Users\User\Desktop\research\multiviewlearning\shiguoxin\bbc_seg14of4.mat');
data=f.data;
label=f.labels;
% addpath('C:\Users\User\Desktop\research\kernelclusteringexp')
para1=[10 ];
para2=[.001  ];
para3=[ .00001  ];

% for i=1:size(data,1)
% dist = max(max(data{i})) - min(min(data{i}));
% m01 = (data{i} - min(min(data{i})))/dist;
% data{i} = 2 * m01 - 1;
% end
B=cell(size(data));
for i=1:length(para1)
   for ii=1:size(data,1)
           B{ii}= inv(data{ii}*data{ii}'+para1(i)*eye(size(data{1},1)));
   end
    for j=1:length(para2)
        for k=1:length(para3)
%             fprintf('params%12.6f%12.6f%12.6f\n',para1,para2,para3)
            [result,YY,Fv,S]=secondnorm22(data,label,B,para1(i),para2(j),para3(k));
            dlmwrite('bbc.txt',[para1(i) para2(j) para3(k) result(1,:) result(2,:) result(3,:)   ],'-append','delimiter','\t','newline','pc');
        end
    end
end
        