% 生成模拟星图

% %% 读取数据文件并进行预处理
% 
% 
% data = readtable('data/hip_table.dat');
% % 除去最后一行
% data = data(1:end-1, :);
% % 取第二列作为编号第四列作为数据
% column2 = data{:, 2};
% column4 = data{:, 4};
% column5 = data{:, 5};
% STAR_LENGTH = length(column2); 
% STAR_MAXNUM = max(column2);
% STAR_MAXMAG = 7.5;
% 
% 
% VectorCatalog = cell(STAR_LENGTH,1);
% SphareCatalog = cell(STAR_LENGTH,1);
% % 使用 cellfun 和 strsplit 分割每个 cell 中的字符串  
% splitStrings = cellfun(@strsplit, column4, 'UniformOutput', false);  
% % 遍历分割后的字符串，将它们转换为数值并整形为 3x2 数组  
% for i = 1:STAR_LENGTH
%     % 将字符串转换为数值  
%     coor = str2double(splitStrings{i});  
%       
%     % 将右升交点从小时转换为度  RAhms
%     longitude = (coor(1) + coor(2)/60 + coor(3)/3600) * (360/24); % 1小时 = 15度  
%     % 将赤纬从度、分、秒转换为度  DEdms
%     latitude = coor(4) + coor(5)/60 + coor(6)/3600;  
%       
% %     % 转换为笛卡尔坐标  
% %     x = cosd(longitude) * cosd(latitude);  
% %     y = sind(longitude) * cosd(latitude);  
% %     z = sind(latitude);  
% %     % 归一化为单位向量  
% %     V = [x; y; z] / sqrt(x^2 + y^2 + z^2); 
%     [x, y, z] = sphare2vector(longitude,latitude);
%     V = [x; y; z];
%     S = [longitude; latitude];
%     
%     % 存储结果  
%     VectorCatalog{column2(i)} = {V,column5(i)};
%     SphareCatalog{column2(i)} = {S,column5(i)};   
%     
% end
% 
% save('data/VectorCatalog.mat','VectorCatalog');
% save('data/SphareCatalog.mat','SphareCatalog');

%% 找到世界坐标系中，相机视角内的恒星
load('data/SphareCatalog.mat', 'SphareCatalog');

STAR_LENGTH = 118218;
STAR_MAXNUM = 120416;
STAR_MAXMAG = 7.5;

phi_Ocam = rand()*360;
theta_Ocam = rand()*180;

fprintf('视野矢量的方位角和极角是（ %d ， %d ）\n', phi_Ocam, theta_Ocam);  
[a, b, c] = sphare2vector(phi_Ocam,90-theta_Ocam);

% 原z轴向量  
Z_axis = [0; 0; 1];  
Z_cam_axis = [a; b; c];  
% 计算旋转轴（叉积）  
rotation_axis = cross(Z_axis, Z_cam_axis);  
rotation_axis = rotation_axis / norm(rotation_axis);

% 计算旋转角度（点积和反正切）  
theta = acos(dot(Z_axis , Z_cam_axis));  

K = [0, -rotation_axis(3), rotation_axis(2);  
     rotation_axis(3), 0, -rotation_axis(1);  
     -rotation_axis(2), rotation_axis(1), 0];  
matrix_R = eye(3)+sin(theta)*K+(1-cos(theta))*(K^2);
matrix_R = inv(matrix_R);


VisibleStar = [];
% 根据编号遍历数据  
for i = 1:STAR_MAXNUM
    if ~isempty(SphareCatalog{i})  
        phi_P = SphareCatalog{i}{1}(1);
        theta_P = 90 - SphareCatalog{i}{1}(2);
        if SphareCatalog{i}{2} < STAR_MAXMAG && ...
        sind(theta_Ocam)*sind(theta_P)*cosd(phi_Ocam-phi_P) + ...
        cosd(theta_Ocam)*cosd(theta_P) > cosd(6)
        % 球面三角中的余弦定理：cos a = cos b * cos c + sin b sin c cos A
           VisibleStar = [VisibleStar;i,phi_P,90 - theta_P,SphareCatalog{i}{2}];
        end
    else  
%         fprintf('编号 %d 的位置为空\n', i);  
    end  
end



%% 根据视角中的星从世界坐标系转到星敏坐标系

senser_f = 4096*1.414/2/tand(6);
matrix_A = [senser_f,0,2048;
            0,senser_f,2048;
            0,0,1];
        
valid_coords = [];  
valid_xyz = [];  

for i = 1:length(VisibleStar)  
    longitude = VisibleStar(i,2);
    latitude = VisibleStar(i,3);
    [Xstar, Ystar, Zstar] = sphare2vector(longitude,latitude);
    XYZimg = matrix_A * matrix_R * [Xstar; Ystar; Zstar];
    Ximg  = XYZimg(1) / XYZimg(3);
    Yimg  = XYZimg(2) / XYZimg(3);
    if Ximg<0 || Ximg>4096 || Yimg<0 || Yimg>4096
        continue
    end
    valid_xyz = [valid_xyz; Xstar, Ystar, Zstar];
    valid_coords = [valid_coords; Ximg, Yimg];  
end  



%% 生成图像

flag = 2;
% 设置高斯分布的参数
sigma = 12; % 标准差，控制高斯分布的扩散程度
amp = 255;  % 像素值的最大值，这里假设为白色 
neighborhood_size = 10;  

% 创建一个 4096x4096 的黑色背景图像  
img = zeros(4096);  
% 遍历 valid_coords 中的每个坐标点  
for i = 1:size(valid_coords, 1)  
    % 获取当前坐标点  
    x = round(valid_coords(i, 1));  
    y = round(valid_coords(i, 2));  
    
    Kc = 200/2.512^VisibleStar(i,4);
    % 生成高斯分布的二维像素值
    [X, Y] = meshgrid(1:4096, 1:4096);
    if flag == 1
        spot = Kc * amp * exp(-((X - x).^2 + (Y - y).^2) / (2 * sigma^2));
    elseif flag == 2
        spot = 255*(((X - x).^2 + (Y - y).^2) < neighborhood_size^2);
    end
    img = img + spot;
end  


% 添加白噪声
whitenoise = randn(4096) * 0.3; % 生成服从正态分布的随机数
img = img + whitenoise;
  
% 显示图像  
imshow(uint8(img));  
  
% 如果需要，保存图像  
imwrite(uint8(img), 'StarMapSim.png');












    