%%
%Metropolis-Hastings �㷨����ÿ�ε����У���������һ���µ�·����Ȼ����㵱ǰ·������·���ľ�����
% ����·���ľ�����̣���ֱ�ӽ��ܣ�������һ���ĸ��ʽ�����·����
%��������뵱ǰ�¶��йء����ŵ������������ӣ��¶Ȳ����½���������·���ĸ���Ҳ���𽥼�С������㷨���ջ�������һ�����ŵĽ⡣
%%
% ��ʼ������
clc;
clear;
load('data2.mat')
C = data2;
n= size(C,1); %TSP ����Ĺ�ģ,��������Ŀ
D = zeros(n); %�����������о���������
for i = 1:n
    for j = 1:n
        D(i,j) = ((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5;
    end
end
cities=C;
distances=D;
T = 100; % ��ʼ�¶�
Tmin = 1e-3; % ��ֹ�¶�
alpha = 0.99; % ����ϵ��
iter = 1000; % ��������
% ������ɳ�ʼ��
cur_path = randperm(n);
cur_dist = get_dist(cur_path, distances);

% ��������
while T > Tmin
    for i = 1:iter
        % �����½�
        new_path = get_new_path(cur_path);
        new_dist = get_dist(new_path, distances);

        % ����������
        delta_E = new_dist - cur_dist;

        % �ж��Ƿ�����½�
        if delta_E < 0 || exp(-delta_E/T) > rand()
            cur_path = new_path;
            cur_dist = new_dist;
        end
    end
    
    % ����
    T = T * alpha;
end

path_str = num2str(cur_path(1));
for i = 2:length(cur_path)
    path_str = strcat(path_str, '->', num2str(cur_path(i)));
end
disp(['����·��Ϊ��', path_str]);

% ����·��ͼ
figure;
hold on;
for i = 1:n-1
    x = [cities(cur_path(i), 1), cities(cur_path(i+1), 1)];
    y = [cities(cur_path(i), 2), cities(cur_path(i+1), 2)];
    plot(x, y, 'b-');
end
x = [cities(cur_path(n), 1), cities(cur_path(1), 1)];
y = [cities(cur_path(n), 2), cities(cur_path(1), 2)];
plot(x, y, 'b-');
scatter(cities(:, 1), cities(:, 2), 'r', 'filled');
title(['Minlength:',num2str(cur_dist)]);
%�����ܵ�·������
function dist = get_dist(path, distances)
    n = length(path);
    dist = 0;
    for i = 1:n-1
        dist = dist + distances(path(i), path(i+1));
    end
    dist = dist + distances(path(n), path(1));
end
%�����µĽ�
function new_path = get_new_path(path)
    n = length(path);
    i = randi(n-1);
    j = randi(n-1);
    if i == j
        j = mod(j, n-1) + 1;
    end
    if i > j
        tmp = i;
        i = j;
        j = tmp;
    end
    new_path = [path(1:i-1), fliplr(path(i:j)), path(j+1:end)];
end
