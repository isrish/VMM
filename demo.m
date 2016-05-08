close all;
clear;
%% Load data
load('demo_data.mat')
[N,dim] = size(data);

obj=fitVMM_CEM(data,10,'Regularize',0,'debg',1,'tol',1e-13);

% Visualize the result;
figh = figure('Position',[100 950 600 600]);
% www.mathworks.com/matlabcentral/fileexchange/29702-generate-maximally-perceptually-distinct-colors
colors = distinguishable_colors(obj.NComponents*2,'w'); 
for k=1:obj.NComponents
    mu = obj.mu(k,:);
    idx = obj.Class==k;
    if dim>=3
        plot3(data(idx,1),data(idx,2),data(idx,3),'.','MarkerSize', 4, 'Color', colors(k,:), 'Marker', '*');
        hold on;
        h = line([0 mu(1)], [0 mu(2)], [0 mu(3)], 'LineWidth', 2, 'Color', colors(k,:));
    elseif dim==2
        plot(data(idx,1),data(idx,2),'.','MarkerSize', 4, 'Color', colors(k,:), 'Marker', '*');
        hold on;
        h = line([0 mu(1)], [0 mu(2)], 'LineWidth', 2, 'Color', colors(k,:));
    end
end
if dim>=3
    %% Overlay a unit sphere
    % Generate the x, y, and z data for the sphere
    r = 1;
    rg = r * ones(20, 20); % radius is 1
    [th, phi] = meshgrid(linspace(0, 2*pi, 20), linspace(-pi, pi, 20));
    [x,y,z] = sph2cart(th, phi, rg);
    x = x + 0;  % center at 0 in x-direction
    y = y + 0;  % center at 0 in y-direction
    z = z + 0;   % center at0 in z-direction
    lightGrey = 0.8*[1 1 1]; % It looks better if the lines are lighter
    surface(x,y,z,'FaceColor', 'none','EdgeColor',lightGrey); hold on;
    box on;
    grid on;
    view([1 1 0.75]); % adjust the viewing angle
    set(gca,'xlim',[-r,r],'ylim',[-r,r],'zlim',[-r,r]);
end

