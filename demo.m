clc;
clear;

%% 加载数据
load 3sources
Y=truelabel{1}';
nCluster=length(unique(Y));

%% 数据预处理
A=[];
for iv = 1:3
    X{iv}=data{iv};
    X{iv} = NormalizeFea(X{iv},0);
    A=[A;X{iv}];
end
A = NormalizeFea(A,0);
stringaaa ={  'spearman', 'correlation', 'cosine'  ,'jaccard','euclidean'  };
for i=1:length(stringaaa)
    B{i}=pdist2( A' , A',stringaaa{i});
    D{i}=B{i};
end


alf= 4.9;
beita= 2;
gamma= 0.3;


%% run our methods
[P1]=MVC3MF(D,Y, alf,beita, gamma );
for rep=1:10
    [nmi(rep),acc(rep),f(rep),purity(rep),idx]=zhixing_Kmeans(P1,nCluster,Y);
end
ACC=mean(acc);
NMI=mean(nmi);
F=mean(f);
Purity=mean(purity);
[ACC,NMI,F,Purity]
