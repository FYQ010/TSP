function out= loc(data)
% 地球的半径
R = 6371;
% 中央经线的经度
lambda0 = 0;
% 经纬度坐标
lon =data(:,1);
lat =data(:,2);
% 将经纬度转换为笛卡尔坐标
out(:,1)= R * deg2rad(lon - lambda0);
out(:,2) = R * log(tan(pi/4 + deg2rad(lat)/2));
end