function [E,F,G,H]=diophantine(A,B,j)
% j=5;
% A=[1 -1.5 0.7];
% B=[1 0.5];
na=length(A);
nb=length(B);
E(1,:)=[1,zeros(1,j-1)];
n=length(conv(A,[1 -1]));
Q=[1,zeros(1,n-1)];
F1=Q-conv(A,[1 -1]);
n=length(F1);
F(1,:)=F1(2:n);
G(1,:)=[B(1),zeros(1,j-1)];
H(1,:)=B(2);
for i=2:j
    F1=F(i-1,1)*[zeros(1,i-1),1];
    n=length(F1);
    E(i,:)=E(i-1,:)+[F(i-1,1)*[zeros(1,i-1),1],zeros(1,j-n)];
    Fl=[F(i-1,:),0]-F(i-1,1)*conv(A,[1 -1]);
    n1=length(Fl);
    F(i,:)=Fl(2:n1);
end
for t=2:j
    T=conv(E(t,:),[B,zeros(1,j-nb)]);
    H(t,:)=T(t+1);
    G(t,:)=[T(1:t),zeros(1,j-length(T(1:t)))];
end