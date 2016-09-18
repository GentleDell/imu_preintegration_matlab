classdef KITTIGT
    
    methods(Static = true)
        function position = LoadTrajactory( fileName )
            %LOADTXT 此处显示有关此函数的摘要
            %   此处显示详细说明
            aa = load(fileName);
            [rows cols] = size(aa);
            % pos = zeros(rows,3);
            position = [aa(:,4),aa(:,8),aa(:,12)];
            % for i = 1:rows
            %    T = zeros(3,4);
            %    T(1,:) = aa(i,1:4);
            %    T(2,:) = aa(i,5:8);
            %    T(3,:) = aa(i,9:12);
            %
            % end
        end
        
        function pose = LoadPose( fileName )
            aa = load(fileName);
            [rows cols] = size(aa);
            pose = zeros(3,4,rows);
            
            for i = 1:rows
                pose(1,:,i) = aa(i,1:4);
                pose(2,:,i) = aa(i,5:8);
                pose(3,:,i) = aa(i,9:12);
                
            end
        end
    end
    
end


