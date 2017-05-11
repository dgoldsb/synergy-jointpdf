% This script scores models for mets by ChEBI IDs included.

all_IDs = union(iFF.metChEBIID, iND.metChEBIID);
all_IDs = union(all_IDs, iIN.metChEBIID);
all_IDs = union(all_IDs, iMM.metChEBIID);
all_IDs = union(all_IDs, Y4.metChEBIID);
all_IDs = union(all_IDs, iAZ.metChEBIID);
all_IDs = union(all_IDs, bs.metChEBIID);
all_IDs = union(all_IDs, Y5.metChEBIID);
all_IDs = union(all_IDs, iTO.metChEBIID);
all_IDs = union(all_IDs, Y6.metChEBIID); % a total of 1141 ids
all_IDs = union(all_IDs, Y7.metChEBIID); % 1141 ids
all_IDs = union(all_IDs, bio.metChEBIID); % 1718 ids


ID_scores = zeros(length(all_IDs), 12);
ID_scores(:, 1) = ismember(all_IDs,iFF.metChEBIID);
ID_scores(:, 2) = ismember(all_IDs,iND.metChEBIID);
ID_scores(:, 3) = ismember(all_IDs,iIN.metChEBIID);
ID_scores(:, 4) = ismember(all_IDs,iMM.metChEBIID);
ID_scores(:, 5) = ismember(all_IDs,Y4.metChEBIID);
ID_scores(:, 6) = ismember(all_IDs,iAZ.metChEBIID);
ID_scores(:, 7) = ismember(all_IDs,bs.metChEBIID);
ID_scores(:, 8) = ismember(all_IDs,Y5.metChEBIID);
ID_scores(:, 9) = ismember(all_IDs,iTO.metChEBIID);
ID_scores(:, 10) = ismember(all_IDs,Y6.metChEBIID);
ID_scores(:, 11) = ismember(all_IDs,Y7.metChEBIID);
ID_scores(:, 12) = ismember(all_IDs,bio.metChEBIID);

ID_sums = sum(ID_scores,2);
normalized_ID_sums = ID_sums/max(ID_sums);

% visualization
labels = {'iFF', 'iND', 'iIN', 'iMM', 'Y4', 'iAZ', 'iMMbs', 'Y5', ...
    'iTO', 'Y6', 'Y7', 'Bio'};

% test = clustergram(ID_scores, 'ColumnLabels', labels);
% 
% Y = pdist(ID_scores');
% [Z,eigvals] = cmdscale(Y);
% 
% 
% % displacement so the text does not overlay the data points
% dx = 1; 
% dy = 1; 
% 
% scatter3(Z(:,1),Z(:,2),Z(:,3));
% text(Z(:,1)+dx, Z(:,2)+dy, Z(:,3), labels);

% visualization, courtesy of James Eddy

met_clustergram = ...
    clustergram(ID_scores, 'ColumnLabels', labels);
permRows = str2num(cell2mat(get(met_clustergram,'RowLabels')));

close all

% Generate hierarchical clustering distance & linkage to make dendrogram
Y = pdist(ID_scores');
Z = linkage(Y);

% Plot dendrogram in upper 20% of the figure
subplot('position',[.05,.8,.9,.2])
[~,~,permCols] = dendrogram(Z); % save column indices to re-order heatmap
xlim([.5 12.5]) % change x-axis limits to line up with heatmap

% Adjust a few other appearance settings
set(gca,'ytick',[],'color','none','visible','off')
box off

% Create custom colormap
yellow = [255,192,0]/255;
scale1 = linspace(0,yellow(1),128);
scale2 = linspace(0,yellow(2),128);
scale3 = linspace(0,yellow(3),128);
yellowmap = [scale1',scale2',scale3'];
colormap(yellowmap)

% Plot heatmap in lower ~70% of the figure; save room for labels
subplot('position',[.05,.1,.9,.7])
imagesc(ID_scores(flipud(permRows),permCols))

% Adjust appearance and add column labels
set(gca,'ytick',[],'xtick',[])
ylabel('Metabolites','fontsize',8,'fontweight','bold')
for i = 1:12
    text(i,size(ID_scores,1)+48,labels(permCols(i)),...
        'fontsize',8,'horizontalalignment','right','rotation',45)
end

% Create a separate figure for MDS scatterplot
figure;
[Z,eigvals] = cmdscale(Y);

% displacement so the text does not overlay the data points
dx = 2; 
dy = 2; 

% Custom colors for markers
blue = [198,217,241]/255;
dark_blue = [31,73,125]/255;

% Plot scaled values, adjust appearance
s3 = scatter3(Z(:,1),Z(:,2),Z(:,3));
set(s3,'markerfacecolor',dark_blue,'markeredgecolor','k')
set(gca,'xticklabel',[],'yticklabel',[],'zticklabel',[])

% Add labels
text(Z(:,1)+dx, Z(:,2)+dy, Z(:,3), labels, 'FontSize', 8);