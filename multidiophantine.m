function [E,F,G]=multidiophantine(a,b,c,N)
% a��A(z^-1)*��=A�������С�=1-z^-1
% b��B(z^-1)=B
% c��1
% N��Ԥ�ⳤ��
na=length(a)-1; % aΪAʱ��na=3
nb=length(b)-1; % bΪBʱ��nb=1
nc=length(c)-1; % cΪ1ʱ��nc=0
% E,F,G�ĳ�ֵ
E=zeros(N);
E(1,1)=1; % E�ĳ�ֵE1=1
F(1,:)=conv(b,E(1,:)); % F1=B*E1
if na>=nc
    G(1,:)=[c(2:nc+1) zeros(1,na-nc)]-a(2:na+1); %��c(nc+2)=c(nc+3)=...=0
else
    G(1,:)=c(2:nc+1)-[a(2:na+1) zeros(1,nc-na)]; %��a(na+2)=a(na+3)=...=0
end
% ��E,G,F
for j=1:N-1
    for i=1:j
        E(j+1,i)=E(j,i);
    end
    E(j+1,j+1)=G(j,1); % ���µ���E
    for i=2:na
        G(j+1,i-1)=G(j,i)-G(j,1)*a(i);
    end
    G(j+1,na)=-G(j,1)*a(na+1); % ���µ���G
    F(j+1,:)=conv(b,E(j+1,:)); % ���µ���F
end