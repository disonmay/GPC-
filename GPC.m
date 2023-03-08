clc;
clear;
N=5; % 预测长度
Nu=3; % 控制长度
L=100; % 仿真步数
alpha=0.7; % 柔化因子
a=[1 -1.5 0.7]; % A(z^-1)=1-1.5z^(-1)+0.7z^(-2)
b=[1 0.5]; % B(z^-1)=1+0.5z^(-1)
guiji=input("请输入给定信号类型：1——方波；2——阶跃；3——正弦；\n");
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
suanfa=input("请给定算法：1——定常GPC；2——自校正GPC；\n");
if suanfa==1
    c = 1;
    d = 1;   % 滞后延迟
    na=length(a)-1;      % A的阶数
    b=[zeros(1,d-1) b];  % d不是1，所以添0
    nb=length(b)-1;      % B的阶数
    aa=conv(a,[1 -1]);    % 等同于多项式的乘积
    naa=na+1;     
    N1=d;  % 最小输出长度
    gamma=1;  % 控制加权矩阵
    uk=zeros(d+nb,1); % 输入初值 uk(i)表示u(k-i)
    duk=zeros(d+nb,1);   % 控制增量初值
    yk=zeros(naa,1);   % 输出初值
    xi=sqrt(0.01)*randn(L,1);  % 白噪声序列
    % 求解多步丢番图方程并构建F1,F2,G
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
    % 开始采集输出值
    for k=1:L
        time(k)=k;
        y(k)=-aa(2:naa+1)*yk+b*duk(1:nb+1)+xi(k);  % 采集当前输出数据
        Yk=[y(k);yk(1:na)];  % 当前与过去的实际输出Y(k)
        dUK=duk(1:nb);        % 构建向量 delta U(k-j) ：过去的控制增量向量
        % 参考轨迹
        yr(k)=y(k);
        for i=1:N
            yr(k+i)=alpha*yr(k+i-1)+(1-alpha)*w(k+1);
        end
        Yr=[yr(k+N1:k+N)]';  % 对未来参考轨迹输出Yr(k)
        % 求控制量
        dU=inv(F1'*F1+gamma)*F1'*(Yr-F2*dUK-G*Yk);  % △U
        du(k)=dU(1);
        u(k)=uk(1)+du(k);
        % 更新数据
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
    title("预测值与给定值");
    subplot(2,1,2);
    plot(time,u);
    xlabel('k');
    ylabel('u(k)');
    title("控制量u");
else
    c = 1;
    d = 1;   % 滞后延迟
    na=length(a)-1;      % A的阶数
    b=[zeros(1,d-1) b];  % d不是1，所以添0
    nb=length(b)-1;      % B的阶数
    naa=na+1;     
    N1=d;  % 最小输出长度
    gamma=1*eye(Nu);  % 控制加权矩阵
    uk=zeros(d+nb,1); % 输入初值 uk(i)表示u(k-i)
    duk=zeros(d+nb,1);   % 控制增量初值
    yk=zeros(na,1);   % 输入初值
    dyk=zeros(na,1);  % 输入增量初值
    xi=sqrt(0.01)*randn(L,1);  % 白噪声序列
    % RLS初值
    thetae_1=0.001*ones(na+nb-d+2,1);  % 不辨识b中添加的0
    P=10^6*eye(na+nb-d+2);
    lambda=1;  % 遗忘因子[0.9 1]
    for k=1:L
        time(k)=k;
        dy(k)=-a(2:na+1)*dyk(1:na)+b*duk(1:nb+1)+xi(k);  
        y(k)=yk(1)+dy(k);  % 采集当前输出数据
        Yk=[y(k);yk(1:na)];  % 当前与过去的实际输出Y(k)
        dUK=duk(1:nb);        % 构建向量 △U(k-j) ：过去的控制增量向量
        % 参考轨迹
        yr(k)=y(k);
        for i=1:N
            yr(k+i)=alpha*yr(k+i-1)+(1-alpha)*w(k+d);
        end
        Yr=[yr(k+N1:k+N)]';  % 对未来参考轨迹输出Yr(k)
        % 递推最小二乘法
        phi=[-dyk(1:na);duk(d:nb+1)];
        K=P*phi/(lambda+phi'*P*phi);
        thetae(:,k)=thetae_1+K*(dy(k)-phi'*thetae_1);
        P=(eye(na+nb-d+2)-K*phi')*P/lambda;
        % 提取辨识参数
        ae=[1 thetae(1:na,k)'];
        be=[zeros(1,d-1) thetae(na+1:na+nb-d+2,k)'];
        aae=conv(ae,[1 -1]);
        % 求解多步丢番图方程并构建F1,F2,G
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
        % 求控制量
        dU=inv(F1'*F1+gamma)*F1'*(Yr-F2*dUK-G*Yk);  % △U
        du(k)=dU(1);
        u(k)=uk(1)+du(k);
        % 更新数据
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
    title("预测值与给定值");
    subplot(2,1,2);
    plot(time,u);
    xlabel('k');
    ylabel('u(k)');
    axis([0 L -4 4]);
    title("控制量u");
    figure(2);
    plot([1:L],thetae);
    xlabel('k');
    ylabel('辨识参数a,b');
    legend('a_1','a_2','b_0','b_1');
    axis([0 L -3 3]);
    title("模型辨识参数值");
end