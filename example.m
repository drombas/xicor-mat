% Example of Chaterjee's xi-correlation.
%
% David Romero-Bascones (dromero@mondragon.edu)
% Biomedical Engineering Department, Mondragon Unibertsitatea, 2022.

close all;clc;clearvars;

n_example = 8;

n = 100;
x = linspace(-3,3,n);

% Linear relationship
Y(1,:) = x + 0.2*randn(1,n);
Y(2,:) = x + randn(1,n);
Y(3,:) = x + 5*randn(1,n);
Y(4,:) = -x + 0.2*randn(1,n);

% Quadratic relationship
Y(5,:) = x.^2 + 0.5*randn(1,n);
Y(6,:) = x.^2 + 2*randn(1,n);

% Oscilatory relationship
Y(7,:) = cos(2*x) + 0.1*randn(1,n);
Y(8,:) = cos(2*x) + 0.5*randn(1,n);

figure('Position',[500 500 1200 600]);
colors = rand(n_example, 3);
for i=1:n_example
    y = Y(i,:);
    
    % Pearson correlation
    r = corrcoef(x,y);
    r = r(1,2);
    
    % Chaterjee's xi-correlation
    xi = xicor(x,y);
    
    % rRound values
    r = round(r,2);
    xi = round(xi,2);
    
    % Visualization
    subplot(2,4,i);
    scatter(x,y,'filled','MarkerFaceColor',colors(i,:));
    xlim([-3.2 3.2]);
    ylim([min(y)-0.6 max(y)+0.6]);

    xticklabels([]);
    yticklabels([]);
    xticks([]);
    yticks([]);
    box on;
    h = gca;
    h.LineWidth=2;    
    title(['r:' num2str(r) '   \xi_{n}:' num2str(xi)]);
    set(gca,'fontname','Lato');
end
