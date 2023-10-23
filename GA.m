clear;clc;
clc; %清屏
load('data2.mat')
C = data2;
N = size(C,1); %TSP 问题的规模,即城市数目
D = zeros(N); %任意两个城市距离间隔矩阵
for i = 1:N
    for j = 1:N
        D(i,j) = ((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5;
    end
end
NP =1200; %种群规模
G = 600; %最大遗传代数
f = zeros(NP,N); %用于存储种群
F = []; %种群更新中间存储
for i = 1:NP
    f(i,:) = randperm(N); %随机生成初始种群
end
R = f(1,:); %存储最优种群
len = zeros(NP,1); %存储路径长度
fitness = zeros(NP,1); %存储归一化适应值
gen = 0;
%主循环
while gen < G
    %路径长度
    for i = 1:NP
        len(i,1) = D(f(i,N),f(i,1));
        for j = 1:(N-1)
            len(i,1) = len(i,1)+D(f(i,j),f(i,j+1));
        end
    end
    maxlen = max(len); %最长路径
    minlen = min(len); %最短路径
    %更新最短路径
    rr = find(len==minlen);
    R = f(rr(1,1),:);
    %计算归一化适应值
    for i = 1:length(len)
        fitness(i,1) = (1-((len(i,1)-minlen)/(maxlen-minlen+0.001)));
    end
    %选择操作
    nn = 0;
    for i = 1:NP
        if fitness(i,1) >= rand
            nn = nn+1;
            F(nn,:) = f(i,:);
        end
    end
    [aa,bb] = size(F);
    while aa < NP
        nnper = randperm(nn);
        A = F(nnper(1),:);
        B = F(nnper(2),:);
        %交叉操作
        if rand<0.7%交叉概率
            W = ceil(N/10); %交叉点个数
            p = unidrnd(N-W+1); %随机选择交叉范围，从 p 到 p+W
            for i = 1:W
                x = find(A==B(1,p+i-1));
                y = find(B==A(1,p+i-1));
                temp = A(1,p+i-1);
                A(1,p+i-1) = B(1,p+i-1);
                B(1,p+i-1) = temp;
                temp = A(1,x);
                A(1,x) = B(1,y);
                B(1,y) = temp;
            end
        end
        %变异操作
        if rand<0.5%变异概率
            p1 = floor(1+N*rand());
            p2 = floor(1+N*rand());
            while p1==p2
                p1 = floor(1+N*rand());
                p2 = floor(1+N*rand());
            end
            tmp = A(p1);
            A(p1) = A(p2);
            A(p2) = tmp;
            tmp = B(p1);
            B(p1) = B(p2);
            B(p2) = tmp;
        end
        F = [F;A;B];
        [aa,bb] = size(F);
    end
    if aa > NP
        F = F(1:NP,:); %保持种群规模为 NP
    end
    f = F; %更新种群
    f(1,:) = R; %保留每代最优个体
    clear F;
    gen = gen+1;
    Rlength(gen) = minlen;
end
figure
for i = 1:N-1
    plot([C(R(i),1),C(R(i+1),1)],[C(R(i),2),C(R(i+1),2)],'b-');
    hold on;
end
plot([C(R(N),1),C(R(1),1)],[C(R(N),2),C(R(1),2)],'r-');
scatter(C(:, 1), C(:, 2), 'r', 'filled');
title(['优化最短距离:',num2str(minlen)]);
figure
plot(Rlength)
xlabel('迭代次数')
ylabel('目标函数值')
title('适应度进化曲线')

[~,name] = xlsread('中国一级行政区坐标.xlsx','中国一级行政区坐标','A1:A34');
R1 = '16->2->12->11->10->1->13->14->8->7->28->15->6->5->32->30->4->31->9->27->34->33->25->21->26->20->3->18->19->22->23->24->29->17';
C1 = strsplit(R1, '->');
path_str = num2str(name{str2double(C1{1})});
for i = 2:length(R)
    path_str = strcat(path_str, '->', num2str(cur_path(i)));
end
disp(['最优路径为：', path_str]);

