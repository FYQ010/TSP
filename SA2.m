%%
%Metropolis-Hastings 算法。在每次迭代中，我们生成一个新的路径，然后计算当前路径和新路径的距离差。如
% 果新路径的距离更短，则直接接受；否则，以一定的概率接受新路径，
%这个概率与当前温度有关。随着迭代次数的增加，温度不断下降，接受新路径的概率也会逐渐减小，因此算法最终会收敛到一个较优的解。
%%
% 初始化参数
clc;
clear;
load('data2.mat')
C = data2;
n= size(C,1); %TSP 问题的规模,即城市数目
D = zeros(n); %任意两个城市距离间隔矩阵
for i = 1:n
    for j = 1:n
        D(i,j) = ((C(i,1)-C(j,1))^2+(C(i,2)-C(j,2))^2)^0.5;
    end
end
cities=C;
distances=D;
T = 100; % 初始温度
Tmin = 1e-3; % 终止温度
alpha = 0.99; % 降温系数
iter = 1000; % 迭代次数
% 随机生成初始解
cur_path = randperm(n);
cur_dist = get_dist(cur_path, distances);

% 迭代计算
while T > Tmin
    for i = 1:iter
        % 生成新解
        new_path = get_new_path(cur_path);
        new_dist = get_dist(new_path, distances);

        % 计算能量差
        delta_E = new_dist - cur_dist;

        % 判断是否接受新解
        if delta_E < 0 || exp(-delta_E/T) > rand()
            cur_path = new_path;
            cur_dist = new_dist;
        end
    end
    
    % 降温
    T = T * alpha;
end

path_str = num2str(cur_path(1));
for i = 2:length(cur_path)
    path_str = strcat(path_str, '->', num2str(cur_path(i)));
end
disp(['最优路径为：', path_str]);

% 绘制路径图
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
%计算总的路径长度
function dist = get_dist(path, distances)
    n = length(path);
    dist = 0;
    for i = 1:n-1
        dist = dist + distances(path(i), path(i+1));
    end
    dist = dist + distances(path(n), path(1));
end
%生成新的解
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
