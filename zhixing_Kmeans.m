function [nmi,acc,f,Purity,idx]=zhixing_Kmeans(P,nCluster,gt)
    N=size(P,1);
    P1=NormalizeFea(P,1);
    MAXiter = 1000;
    REPlic = 20;
   [idx,Ctrs,SumD,D]=kmeans(P1,nCluster,'maxiter',MAXiter,'replicates',REPlic,'EmptyAction','singleton');
 %[idx,Ctrs,SumD,D]=kmeans(P,nCluster,'emptyaction', 'singleton','replicates', 100, 'display', 'off');
    [A nmi avgent] = compute_nmi(gt,idx);
  
   acc = Accuracy(idx,double(gt));
    [f,p,r] = compute_f(gt,idx);
    
    predLidx = unique(idx); 
    pred_classnum = length(predLidx);
    correnum = 0;
    for ci = 1:pred_classnum
        incluster = gt(find(idx== predLidx(ci)));
        inclunub = hist(incluster, 1:max(incluster)); 
        if isempty(inclunub) inclunub=0;end;
        correnum = correnum + max(inclunub);
    end;
    Purity = correnum/length(idx);


end