function [result,YY,FV,ZV]=secondnorm22(X,s,B,alpha,beta,gamma)
% s is the true class label.
%
[v1,v2]=size(X);
FV=cell(v1,v2);
ZV=cell(v1,v2);
Wv=zeros(v1,1);
LV=cell(v1,1);
c=length(unique(s));
n=size(X{1},1);

for num = 1:v1
    Fv=randn(n,c);
    Fv= orth(Fv);
    FV{num}=Fv;
    Zv=eye(n);
    ZV{num}=Zv;
end
YY = rand(n,c);
YY = orth(YY);

for i=1:200

    Yold=YY;
    for num = 1:v1
        Zv=ZV{num};
        xv=X{num};
        xv=xv';
        Fv=FV{num};
         Wv(num)=1/(2*norm(YY*YY'-Fv*Fv','fro'));
        parfor ij=1:n
               d=distance(Fv,n,ij);
            Zv(:,ij)=B{num}*(xv'*xv(:,ij) - beta/4*d'); 
        end
        Zv(find(Zv<0))=0;
        Zv=(Zv+Zv')/2;
        D = diag(sum(Zv));
        Lv = D-Zv;
        LL=beta*Lv-2*gamma*Wv(num)*YY*YY'+gamma*Wv(num)*eye(n);
        [Fv, temp, ev]=eig1(LL, c, 0); 
        FV{num}=Fv;
        ZV{num}=Zv;
        LV{num}=Lv;        
    end
    
    A=zeros(n);
     for num=1:v1
        Fv=FV{num};
        A=Wv(num)*(eye(n)-2*Fv*Fv')+A;   
     end
   [YY, temp, ev]=eig1(A, c, 0);
             
        if i>5 &&((norm(YY-Yold)/norm(Yold))<1e-3)
            break
        end 
end
for ij=1:10
actual_ids= kmeans(YY, c, 'emptyaction', 'singleton', 'replicates', 1, 'display', 'off');
[res(ij,:)] = ClusteringMeasure( actual_ids,s);
end
result(1,1)=mean(res(:,1));result(1,2)=std(res(:,1));
result(2,1)=mean(res(:,2));result(2,2)=std(res(:,2));
result(3,1)=mean(res(:,3));result(3,2)=std(res(:,3));

end
function [all]=distance(F,n,ij);
  for ji=1:n
            all(ji)=(norm(F(ij,:)-F(ji,:)))^2;
  end
end

