function plotFlow3(W, H, D, u, v, w , rSize, scale)

figure('position', [500 500 1400 900]);

% [x, y, z] = meshgrid(1:W/rSize(1),1:H/rSize(2),1:D/rSize(3));
[x, y, z] = meshgrid(W/rSize(1):-1:1,H/rSize(2):-1:1,D/rSize(3):-1:1);

% quiver3(x, y, z, resize(u, size(x)), resize(v, size(y)), ...
%     resize(w, size(z)), scale, 'color', 'r', 'linewidth', 2);

% u_tmp = resize(u, size(x));
% v_tmp = resize(v, size(x));
% w_tmp = resize(w, size(x));
% 
% color = zeros( [size(u_tmp) 3] );
% for i=1:size(u_tmp,3)
%     color(:,:,i,:) = reshape( hsv2rgb( cat(3,360*sqrt( u_tmp(:,:,i).^2 + v_tmp(:,:,i).^2 + w_tmp(:,:,i).^2 ),...
%         0.5*ones( size( u_tmp(:,:,i)) ), 0.5*ones( size( u_tmp(:,:,i)) ) ) ), 5, 10, 1, 3);
% end
% 
% quiver3(x, y, z, u_tmp, v_tmp, w_tmp, scale, 'color', ...
%     color, 'linewidth', 2);
% set(gca,'YDir','reverse');

%%

% u_tmp = resize(u, size(x));
% v_tmp = resize(v, size(x));
% w_tmp = resize(w, size(x));

q = quiver3(x, y, z, resize(u, size(x)), resize(v, size(y)), ...
    resize(w, size(z)), scale, 'linewidth', 2);
xlabel('z');
ylabel('y');
zlabel('x');

mags = sqrt(sum(cat(2, q.UData(:), q.VData(:), ...
            reshape(q.WData, numel(q.UData), [])).^2, 2));

%// Get the current colormap
currentColormap = colormap(jet);

%// Now determine the color to make each arrow using a colormap
[~, ~, ind] = histcounts(mags, size(currentColormap, 1));

%// Now map this to a colormap to get RGB
cmap = uint8(ind2rgb(ind(:), currentColormap) * 255);
cmap(:,:,4) = 255;
cmap = permute(repmat(cmap, [1 3 1]), [2 1 3]);

%// We repeat each color 3 times (using 1:3 below) because each arrow has 3 vertices
set(q.Head, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:3,:,:), [], 4).');   %'

%// We repeat each color 2 times (using 1:2 below) because each tail has 2 vertices

set(q.Tail, ...
    'ColorBinding', 'interpolated', ...
    'ColorData', reshape(cmap(1:2,:,:), [], 4).');

colorbar
set(gca, 'CLim', [min(mags), max(mags)]);
% set(gca, 'CLim', [min(mags)*10, max(mags)*10]);

h=colorbar;
set(h, 'Position', [.92 .3 .03 .65] );

set(gca,'xtick',[0 2.5 5]);
set(gca,'ytick',[0 1.5 3.0 4.5 6.0]);
set(gca,'ztick',[0 2.5 5]);

set(gca,'xticklabel',[0 50 100]);
set(gca,'yticklabel',[0 45 90 135 180]);
% set(gca,'xticklabel',[100 50 0]);
% set(gca,'yticklabel',[180 135 90 45 0]);
set(gca,'zticklabel',[0 50 100]);

%% real flow
% set(gca,'xtick',[0 3 6 9 12]);
% set(gca,'ytick',[0 1.5 3.0 4.5 6.0]);
% set(gca,'ztick',[0 1.5 3]);
% set(gca,'xticklabel',[0 100 200 300 400]);
% set(gca,'yticklabel',[0 45 90 135 180]);
% set(gca,'zticklabel',[0 40 80]);

end
