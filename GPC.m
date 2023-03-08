clc;
clear;
N=5; % Ԥ�ⳤ��
Nu=3; % ���Ƴ���
L=100; % ���沽��
alpha=0.7; % �ữ����
a=[1 -1.5 0.7]; % A(z^-1)=1-1.5z^(-1)+0.7z^(-2)
b=[1 0.5]; % B(z^-1)=1+0.5z^(-1)
guiji=input("����������ź����ͣ�1����������2������Ծ��3�������ң�\n");
switch guiji
    case 1
       w=10*[ones(L/4,1);-ones(L/4,1);ones(L/4,1);-ones(L/4+1,1)];  
    case 2
       w=10*[ones(L/4,1);ones(L/4,1);ones(L/4,1);ones(L/4+1,1)];   
    case 3
        for i=1:L+1                                               
            p(i)=10*sin(2*pi*i/20);
            w=p;
        end
end
suanfa=input("������㷨��1��������GPC��2������У��GPC��\n");
if suanfa==1
    c = 1;
    d = 1;   % �ͺ��ӳ�
    na=length(a)-1;      % A�Ľ���
    b=[zeros(1,d-1) b];  % d����1��������0
    nb=length(b)-1;      % B�Ľ���
    aa=conv(a,[1 -1]);    % ��ͬ�ڶ���ʽ�ĳ˻�
    naa=na+1;     
    N1=d;  % ��С�������
    gamma=1;  % ���Ƽ�Ȩ����
    uk=zeros(d+nb,1); % �����ֵ uk(i)��ʾu(k-i)
    duk=zeros(d+nb,1);   % ����������ֵ
    yk=zeros(naa,1);   % �����ֵ
    xi=sqrt(0.01)*randn(L,1);  % ����������
    % ���ಽ����ͼ���̲�����F1,F2,G
    [E,F,G]=multidiophantine(aa,b,c,N);
    G=G(N1:N,:);            % F
    F1=zeros(N-N1+1,Nu);    % G
    F2=zeros(N-N1+1,nb);    % H
    for i=1:N-N1+1
        for j=1:min(i,Nu)
        F1(i,j)=F(i+N1-1,i+N1-1-j+1);
        end
        for j=1:nb
        F2(i,j)=F(i+N1-1,i+N1-1+j);
        end
    end
    % ��ʼ�ɼ����ֵ
    for k=1:L
        time(k)=k;
        y(k)=-aa(2:naa+1)*yk+b*duk(1:nb+1)+xi(k);  % �ɼ���ǰ�������
        Yk=[y(k);yk(1:na)];  % ��ǰ���ȥ��ʵ�����Y(k)
        dUK=duk(1:nb);        % �������� delta U(k-j) ����ȥ�Ŀ�����������
        % �ο��켣
        yr(k)=y(k);
        for i=1:N
            yr(k+i)=alpha*yr(k+i-1)+(1-alpha)*w(k+1);
        end
        Yr=[yr(k+N1:k+N)]';  % ��δ���ο��켣���Yr(k)
        % �������
        dU=inv(F1'*F1+gamma)*F1'*(Yr-F2*dUK-G*Yk);  % ��U
        du(k)=dU(1);
        u(k)=uk(1)+du(k);
        % ��������
        for i=1+nb:-1:2
            uk(i)=uk(i-1);
            duk(i)=duk(i-1);
        end
        uk(1)=u(k);
        duk(1)=du(k);
        for i=naa:-1:2
            yk(i)=yk(i-1);
        end
        yk(1)=y(k);
    end
    subplot(2,1,1);
    plot(time,w(1:L),'r:',time,y);
    xlabel('k');
    ylabel('w(k) y(k)');
    legend('w(k)','y(k)');
    title("Ԥ��ֵ�����ֵ");
    subplot(2,1,2);
    plot(time,u);
    xlabel('k');
    ylabel('u(k)');
    title("������u");
else
    c = 1;
    d = 1;   % �ͺ��ӳ�
    na=length(a)-1;      % A�Ľ���
    b=[zeros(1,d-1) b];  % d����1��������0
    nb=length(b)-1;      % B�Ľ���
    naa=na+1;     
    N1=d;  % ��С�������
    gamma=1*eye(Nu);  % ���Ƽ�Ȩ����
    uk=zeros(d+nb,1); % �����ֵ uk(i)��ʾu(k-i)
    duk=zeros(d+nb,1);   % ����������ֵ
    yk=zeros(na,1);   % �����ֵ
    dyk=zeros(na,1);  % ����������ֵ
    xi=sqrt(0.01)*randn(L,1);  % ����������
    % RLS��ֵ
    thetae_1=0.001*ones(na+nb-d+2,1);  % ����ʶb����ӵ�0
    P=10^6*eye(na+nb-d+2);
    lambda=1;  % ��������[0.9 1]
    for k=1:L
        time(k)=k;
        dy(k)=-a(2:na+1)*dyk(1:na)+b*duk(1:nb+1)+xi(k);  
        y(k)=yk(1)+dy(k);  % �ɼ���ǰ�������
        Yk=[y(k);yk(1:na)];  % ��ǰ���ȥ��ʵ�����Y(k)
        dUK=duk(1:nb);        % �������� ��U(k-j) ����ȥ�Ŀ�����������
        % �ο��켣
        yr(k)=y(k);
        for i=1:N
            yr(k+i)=alpha*yr(k+i-1)+(1-alpha)*w(k+d);
        end
        Yr=[yr(k+N1:k+N)]';  % ��δ���ο��켣���Yr(k)
        % ������С���˷�
        phi=[-dyk(1:na);duk(d:nb+1)];
        K=P*phi/(lambda+phi'*P*phi);
        thetae(:,k)=thetae_1+K*(dy(k)-phi'*thetae_1);
        P=(eye(na+nb-d+2)-K*phi')*P/lambda;
        % ��ȡ��ʶ����
        ae=[1 thetae(1:na,k)'];
        be=[zeros(1,d-1) thetae(na+1:na+nb-d+2,k)'];
        aae=conv(ae,[1 -1]);
        % ���ಽ����ͼ���̲�����F1,F2,G
        [E,F,G]=multidiophantine(aae,be,c,N);
        G=G(N1:N,:);            % F
        F1=zeros(N-N1+1,Nu);    % G
        F2=zeros(N-N1+1,nb);    % H
        for i=1:N-N1+1
            for j=1:min(i,Nu)
                F1(i,j)=F(i+N1-1,i+N1-1-j+1);
            end
            for j=1:nb
                F2(i,j)=F(i+N1-1,i+N1-1+j);
            end
        end
        % �������
        dU=inv(F1'*F1+gamma)*F1'*(Yr-F2*dUK-G*Yk);  % ��U
        du(k)=dU(1);
        u(k)=uk(1)+du(k);
        % ��������
        thetae_1=thetae(:,k);
        for i=1+nb:-1:2
            uk(i)=uk(i-1);
            duk(i)=duk(i-1);
        end
        uk(1)=u(k);
        duk(1)=du(k);
        for i=na:-1:2
            yk(i)=yk(i-1);
            dyk(i)=dyk(i-1);
        end
        yk(1)=y(k);
        dyk(1)=dy(k);
    end
    figure(1);
    subplot(2,1,1);
    plot(time,w(1:L),'r:',time,y);
    xlabel('k');
    ylabel('w(k) y(k)');
    legend('w(k)','y(k)');
    title("Ԥ��ֵ�����ֵ");
    subplot(2,1,2);
    plot(time,u);
    xlabel('k');
    ylabel('u(k)');
    axis([0 L -4 4]);
    title("������u");
    figure(2);
    plot([1:L],thetae);
    xlabel('k');
    ylabel('��ʶ����a,b');
    legend('a_1','a_2','b_0','b_1');
    axis([0 L -3 3]);
    title("ģ�ͱ�ʶ����ֵ");
end