function [P,S,Dx,errDx,errS,errP,Allvalue]=MVC3MF(D,gt, alf,beita,gamma,s )
% æ‡¿Î£∫spearman  Euclidean seuclidean correlation cosine Jaccard

m=size(D,2);
n=size(D{1},2);
c=length(unique(gt));
Convalue=[];
P=zeros(n,c);
Dx=zeros(n);
Df=zeros(n);
S=zeros(n);
tempL=zeros(n);
I=2*ones(n,n);
for i=1:m
    w(i)=1/m;
    
end
iter=0;
iterMax=50;
Allvalue=0;
lastDx=zeros(n);
lastS=zeros(n);
lastP=zeros(n,c);
while iter<iterMax
    %  fprintf("ieter=%d  ",iter);
    iter=iter+1 ;
    tempdv=zeros(n);
    for i=1:m
        tempdv=tempdv+2*gamma*w(i)*D{i};
    end
    Dx=tempdv./(I*gamma+2*S);
    
    errDx(iter)=norm( Dx-lastDx,inf)^2;
    lastDx=Dx;
    tempz=Dx.^2+beita*Df;
    S=SimplexProj( -0.5*tempz/alf );
    S= (S'+S)/2 ;
    errS(iter)=norm( S-lastS,inf);
    lastS=S;
    Diag=diag(sum(S));
    L=Diag-S;
    
    %     [uN,sN,vN] = svd(L);
    %     P=uN(:,n-c+1:n);
    [P,~,~]=eig1(L,c,0);
    errP(iter)=norm( P-lastP,inf);
    lastP=P;
    Df=L2_distance_1(P',P');
    
    
    for i=1:m
        w(i)=1/(norm(D{i}-Dx,'fro') );
    end
    sumw=sum(w);
    for i=1:m
        w(i)=w(i)/sumw;
    end
    rep=iter;
    acc(rep)=0;
    vaule=0;
    for i=1:m
        vaule= vaule+gamma*w(i)*norm(D{i}-Dx,'fro').^2;
    end
    vaule=vaule+alf*norm(S,'fro').^2+2*gamma*trace(P'*L*P);
    temp=Dx.^2;
    tempp=temp.*S;
    vaule=vaule+sum(tempp(:));
    Allvalue(iter)=vaule;
    
end
fprintf("\n");
end
