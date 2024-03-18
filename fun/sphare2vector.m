function [x,y,z] = sphare2vector(longitude,latitude)
    % 转换为笛卡尔坐标  
    x = cosd(longitude) * cosd(latitude);  
    y = sind(longitude) * cosd(latitude);  
    z = sind(latitude);  
    
    % 归一化为单位向量  
    magnitude = sqrt(x^2 + y^2 + z^2);  
    V = [x; y; z] / magnitude; % 直接得到归一化后的 3x1 列向量
end

