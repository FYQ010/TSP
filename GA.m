clear;clc;
clc; %����
load('data2.mat')
C = data2;
N = size(C,1); %TSP ����Ĺ�ģ,��������Ŀ
D = zeros(N); %�����������о���������
for i = 1:N
    for j = 1:N
        D(i,j) = ((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5;
    end
end
NP =1200; %��Ⱥ��ģ
G = 600; %����Ŵ�����
f = zeros(NP,N); %���ڴ洢��Ⱥ
F = []; %��Ⱥ�����м�洢
for i = 1:NP
    f(i,:) = randperm(N); %������ɳ�ʼ��Ⱥ
end
R = f(1,:); %�洢������Ⱥ
len = zeros(NP,1); %�洢·������
fitness = zeros(NP,1); %�洢��һ����Ӧֵ
gen = 0;
%��ѭ��
while gen < G
    %·������
    for i = 1:NP
        len(i,1) = D(f(i,N),f(i,1));
        for j = 1:(N-1)
            len(i,1) = len(i,1)+D(f(i,j),f(i,j+1));
        end
    end
    maxlen = max(len); %�·��
    minlen = min(len); %���·��
    %�������·��
    rr = find(len==minlen);
    R = f(rr(1,1),:);
    %�����һ����Ӧֵ
    for i = 1:length(len)
        fitness(i,1) = (1-((len(i,1)-minlen)/(maxlen-minlen+0.001)));
    end
    %ѡ�����
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
        %�������
        if rand<0.7%�������
            W = ceil(N/10); %��������
            p = unidrnd(N-W+1); %���ѡ�񽻲淶Χ���� p �� p+W
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
        %�������
        if rand<0.5%�������
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
        F = F(1:NP,:); %������Ⱥ��ģΪ NP
    end
    f = F; %������Ⱥ
    f(1,:) = R; %����ÿ�����Ÿ���
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
title(['�Ż���̾���:',num2str(minlen)]);
figure
plot(Rlength)
xlabel('��������')
ylabel('Ŀ�꺯��ֵ')
title('��Ӧ�Ƚ�������')

[~,name] = xlsread('�й�һ������������.xlsx','�й�һ������������','A1:A34');
R1 = '16->2->12->11->10->1->13->14->8->7->28->15->6->5->32->30->4->31->9->27->34->33->25->21->26->20->3->18->19->22->23->24->29->17';
C1 = strsplit(R1, '->');
path_str = num2str(name{str2double(C1{1})});
for i = 2:length(R)
    path_str = strcat(path_str, '->', num2str(cur_path(i)));
end
disp(['����·��Ϊ��', path_str]);

