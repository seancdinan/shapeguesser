% Determine the probable shape of an object, and the coordinates of
% its coordinates
% (c) Sean Dinan 2014

% -----------------------------------------------------------------
% Create points
clear; clc;
point_amt1 = 1000; %number of points per region

% Region 1: xmin/max = 0/10 , ymin/max = 0,3
xmin_reg1 = 0;
xmax_reg1 = 10;

ymin_reg1 = 0;
ymax_reg1 = 3;

x_reg1 = xmin_reg1 + (xmax_reg1 - xmin_reg1)*rand(point_amt1,1);
y_reg1 = ymin_reg1 + (ymax_reg1 - ymin_reg1)*rand(point_amt1,1);

xy_reg1 = [x_reg1 , y_reg1];

% Region 2: xmin/max = 0/3 , ymin/max = 3,7
xmin_reg2 = 0;
xmax_reg2 = 3;

ymin_reg2 = 3;
ymax_reg2 = 7;

x_reg2 = xmin_reg2 + (xmax_reg2 - xmin_reg2)*rand(point_amt1,1);
y_reg2 = ymin_reg2 + (ymax_reg2 - ymin_reg2)*rand(point_amt1,1);

xy_reg2 = [x_reg2 , y_reg2];

% Region 3: xmin/max = 0/10 , ymin/max = 7,10
xmin_reg3 = 0;
xmax_reg3 = 10;

ymin_reg3 = 7;
ymax_reg3 = 10;

x_reg3 = xmin_reg3 + (xmax_reg3 - xmin_reg3)*rand(point_amt1,1);
y_reg3 = ymin_reg3 + (ymax_reg3 - ymin_reg3)*rand(point_amt1,1);

xy_reg3 = [x_reg3 , y_reg3];

% Region 4: xmin/max = 7/10 , ymin/max = 3,7
xmin_reg4 = 7;
xmax_reg4 = 10;

ymin_reg4 = 3;
ymax_reg4 = 7;

x_reg4 = xmin_reg4 + (xmax_reg4 - xmin_reg4)*rand(point_amt1,1);
y_reg4 = ymin_reg4 + (ymax_reg4 - ymin_reg4)*rand(point_amt1,1);

xy_reg4 = [x_reg4 , y_reg4];
 
% Combine all regions together
 x_all = [x_reg1;x_reg2;x_reg3;x_reg4];
 y_all = [y_reg1;y_reg2;y_reg3;y_reg4];
 
 xy_all = [x_all, y_all];

% -----------------------------------------------------------------
% Construct the sector grid from mins and maxs. 
% Divide the points into a 16 x 16 grid

n = 2; % Number of sectors in the x or y direction

[sector_max1] = max(xy_all);
[sector_min1] = min(xy_all);

sector_h1 = (sector_max1(2) - sector_min1(2))/n;
sector_w1 = (sector_max1(1) - sector_min1(1))/n;


sector_xy1 = zeros(((n*2)+1),3);

for b = 0:n;
    sector_xy1((b+1),1) = b; % Counter
    sector_xy1((b+1),2) = sector_min1(1) + b*sector_w1; % X
    sector_xy1((b+1),3) = sector_min1(2) + b*sector_h1; % Y
end

fits_trix = zeros(20,2,(n^2)); %used to be (20,2,256)
z = 1;

for i = 1:n;
    for j = 1:n;
        [rws, cls, vls] = find(...
            sector_xy1(i,2) < xy_all(:,1) & ...
            xy_all(:,1) < sector_xy1(i+1,2) & ...
            sector_xy1(j,3) < xy_all(:,2) & ...
            xy_all(:,2) < sector_xy1(j+1,3));
        for k = 1:length(rws);
        fits_trix(k,1,z) = xy_all(rws(k),1);
        fits_trix(k,2,z) = xy_all(rws(k),2);
        end
        z = z+1;
    end
    
% disp(['AFTER J LOOP:', ' i = ',num2str(i),', j = ',num2str(j),', z = ',num2str(z)])
end
z = z - 1; % Compensate for extra addition at end of for loop.

mean1 = zeros(z,2);

for p = 1:z;
    point_amt = 0;
    for q = 1:length(fits_trix(:,:,p));
        if fits_trix(q,1,p) ~= 0 && fits_trix(q,2,p) ~=0;
            point_amt = point_amt + 1;
        end
    end
    if point_amt > 0;
        match_points = fits_trix(1:point_amt,:,p);
        xmean = sum(match_points(:,1)) / point_amt;
        ymean = sum(match_points(:,2)) / point_amt;
        mean1(p,1) = xmean;
        mean1(p,2) = ymean;
    else
        mean1(p,1) = 0;
        mean1(p,2) = 0;
    end
end
mean1 = mean1(any(mean1,2),:);

% Reorganize things to go in order.
ord_mean1 = mean1;
dist = zeros(length(mean1),1);

for s = 2:length(mean1);
    for u = s:length(mean1);
        dist(u,1) = sqrt((ord_mean1(u,1) - ord_mean1((s-1),1))^2 + ...
                         (ord_mean1(u,2) - ord_mean1((s-1),2))^2);
    end
    [minVal,minLoc] = min(dist(s:length(dist)));
    minLoc = minLoc + (s-1);
    holder = ord_mean1(minLoc,:);
    for v = minLoc:-1:s;
        ord_mean1(v,:) = ord_mean1((v-1),:);
    end
    ord_mean1(s,:) = holder;
%     disp(num2str(ord_mean1))
%     disp(' ')
%     disp(' ')
%     pause
end
% disp(['mean1 length: ',num2str(length(mean1))])
% disp(['ord_mean1 length: ',num2str(length(ord_mean1))])
% disp('ord_mean1:')
% disp(num2str(ord_mean1))

plot(xy_all(:,1), xy_all(:,2),'co',...
     ord_mean1(:,1),ord_mean1(:,2),'r-*','LineWidth',1)
    




