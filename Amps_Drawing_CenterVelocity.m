% 清空工作区和图形窗口
clear; close all; clc;

% 加载数据
data1 = load('Eo12.mat');
data2 = load('job2.mat');

% 提取数据并确保为列向量
Time1 = data1.Time(:);  
CenterVelocity1 = data1.CenterVelocity(:);

Time2 = data2.Time(:);  
CenterVelocity2 = data2.CenterVelocity(:);


% 去除零值项（找到最后一个非零值的位置并截断）
idx1 = find(CenterVelocity1 ~= 0, 1, 'last');
Time1 = Time1(1:idx1);
CenterVelocity1 = CenterVelocity1(1:idx1);

idx2 = find(CenterVelocity2 ~= 0, 1, 'last');
Time2 = Time2(1:idx2);
CenterVelocity2 = CenterVelocity2(1:idx2);


% 删除每组数据的最后5项
n = 35;  % 要删除的最后几项
Time1 = Time1(1:end-n-180);
CenterVelocity1 = CenterVelocity1(1:end-n-180);

Time2 = Time2(1:end-n);
CenterVelocity2 = CenterVelocity2(1:end-n);


% 定义颜色和线型
colors = [0.8, 0.2, 0.2;   % 红色
          0.2, 0.2, 0.8;   % 蓝色
          0.2, 0.8, 0.2];  % 绿色
lineStyles = {'-', '--', ':'}; % 不同线型

% 创建图形并绘制所有数据集
figure;
hold on;

% 绘制每组数据
plot(Time1, CenterVelocity1, 'Color', colors(1, :), ...
    'LineWidth', 2, 'LineStyle', lineStyles{1});
plot(Time2, CenterVelocity2, 'Color', colors(2, :), ...
    'LineWidth', 2, 'LineStyle', lineStyles{2});


% 添加标签和标题
xlabel('Non-dimensional Time');
ylabel('Center Velocity');
title('Center Velocity (Y)');

% 添加图例
legend({'Eo = 12', 'Eo = 24'}, 'Location', 'Best');

hold off;
