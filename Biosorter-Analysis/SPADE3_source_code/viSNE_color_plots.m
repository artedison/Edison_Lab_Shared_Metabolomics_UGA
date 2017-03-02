function viSNE_color_plots(cell_positions, data, marker_names, local_density)

[F,XI,hx] = ksdensity(cell_positions(1,:)); 
[F,YI,hy] = ksdensity(cell_positions(2,:)); 
grid_points = [repmat(XI(:)',1,length(YI));reshape(repmat(YI(:)',length(XI),1),1,length(XI)*length(YI))];
window_x = abs(XI(2)-XI(1))*1.1;
window_y = abs(YI(2)-YI(1))*1.1;
freq = zeros(2,length(XI)*length(YI))+NaN;
expr = zeros(length(marker_names),length(XI)*length(YI))+NaN;
for i=1:size(grid_points,2)
    freq(1,i) = sum(exp(-(cell_positions(1,:) - grid_points(1,i)).^2/(2*hx^2)).*exp(-(cell_positions(2,:) - grid_points(2,i)).^2/(2*hy^2))/(sqrt(2*pi*hx^2*2*pi*hy)));
    freq(2,i) = sum((exp(-(cell_positions(1,:) - grid_points(1,i)).^2/(2*hx^2)).*exp(-(cell_positions(2,:) - grid_points(2,i)).^2/(2*hy^2))/(sqrt(2*pi*hx^2*2*pi*hy))).*local_density);
    data_tmp = data(:,abs(cell_positions(1,:) - grid_points(1,i))<=window_x & abs(cell_positions(2,:) - grid_points(2,i))<=window_y);
    if ~isempty(data_tmp)
        expr(:,i) = mean(data_tmp,2);
    end
end
freq(1,:) = freq(1,:)/size(cell_positions,2);
freq(2,:) = freq(2,:)/sum(local_density);


axis_limits = [[min(cell_positions(1,:)), max(cell_positions(1,:))] + [-1,1]*0.05*(max(cell_positions(1,:))-min(cell_positions(1,:))) , [min(cell_positions(2,:)), max(cell_positions(2,:))] + [-1,1]*0.05*(max(cell_positions(2,:))-min(cell_positions(2,:))) ];

% contour plot of downsampled data
h = figure; 

figure(h); subplot(1,2,1); hold off
Z = reshape(freq(1,:),length(XI),length(YI))';
contour(XI,YI,Z,[5:15:60,70:10:90,93:3:99]/100*max(max(Z)));
% hold on; plot(node_positions(1,:),node_positions(2,:),'.b'); drawnow;
% for i=1:size(edges,1), line(node_positions(1,edges(i,:)),node_positions(2,edges(i,:))); end
axis(axis_limits);


% contour plot of original data
figure(h); subplot(1,2,2); hold off
Z = reshape(freq(2,:),length(XI),length(YI))';
contour(XI,YI,Z,[5:15:60,70:10:90,93:3:99]/100*max(max(Z)));
% hold on; plot(node_positions(1,:),node_positions(2,:),'.b'); drawnow;
% for i=1:size(edges,1), line(node_positions(1,edges(i,:)),node_positions(2,edges(i,:))); end
axis(axis_limits);

% contour plot of expr
h = figure;
for i=1:size(expr,1)
    figure(h); 
    subplot_with_k_panels(i,size(expr,1));
    Z = reshape(expr(i,:),length(XI),length(YI))';
    Z(isnan(Z)) = min(expr(i,:)) - (max(expr(i,:))-min(expr(i,:)))*0.01;
    contourf(XI,YI,Z,[1:15:60,70:10:90,95, 99]/100 * (max(expr(i,:))-min(expr(i,:)))+min(expr(i,:)) );
%     contourf(XI,YI,Z,[1, 35, 70, 80, 90, 95, 99]/100 * (max(expr(i,:))-min(expr(i,:)))+min(expr(i,:)) );
    title(marker_names{i});
    axis(axis_limits);
end




function subplot_with_k_panels(i,k)
if k<=3
    subplot(1,k,i);
elseif k<=8
    subplot(2,ceil(k/2),i);
elseif k<=15
    subplot(3,ceil(k/3),i);
else
    subplot(4,ceil(k/4),i);
end
