% This script scores models for genes included.

all_genes = union(iFF.genes, iND.genes);
all_genes = union(all_genes, iIN.genes);
all_genes = union(all_genes, iMM.genes);
all_genes = union(all_genes, Y4.genes);
all_genes = union(all_genes, iAZ.genes);
all_genes = union(all_genes, bs.genes);
all_genes = union(all_genes, Y5.genes);
all_genes = union(all_genes, iTO.genes);
all_genes = union(all_genes, Y6.genes); % a total of 1053 genes
all_genes = union(all_genes, Y7.genes); % 1064 genes

% modify bio model gene annotation
bio.genes  = regexprep(bio.genes,'_i','');
bio.genes  = regexprep(bio.genes,'_\d','');
bio.genes(~cellfun('isempty', (strfind(bio.genes,'MetaCyc')))) = {' '};

all_genes = union(all_genes, bio.genes); % 1187 genes

gene_scores = zeros(length(all_genes), 12);
gene_scores(:, 1) = ismember(all_genes,iFF.genes);
gene_scores(:, 2) = ismember(all_genes,iND.genes);
gene_scores(:, 3) = ismember(all_genes,iIN.genes);
gene_scores(:, 4) = ismember(all_genes,iMM.genes);
gene_scores(:, 5) = ismember(all_genes,Y4.genes);
gene_scores(:, 6) = ismember(all_genes,iAZ.genes);
gene_scores(:, 7) = ismember(all_genes,bs.genes);
gene_scores(:, 8) = ismember(all_genes,Y5.genes);
gene_scores(:, 9) = ismember(all_genes,iTO.genes);
gene_scores(:, 10) = ismember(all_genes,Y6.genes);
gene_scores(:, 11) = ismember(all_genes,Y7.genes);
gene_scores(:, 12) = ismember(all_genes,bio.genes);

gene_sums = sum(gene_scores,2);
normalized_gene_sums = gene_sums/max(gene_sums);

% visualization
labels = {'iFF', 'iND', 'iIN', 'iMM', 'Y4', 'iAZ', 'iMMbs', 'Y5', ...
    'iTO', 'Y6', 'Y7', 'Bio'};

% test = clustergram(gene_scores, 'ColumnLabels', labels, 'FontSize', 18);
% Y = pdist(gene_scores');
% [Z,eigvals] = cmdscale(Y);
% 
% 
% % displacement so the text does not overlay the data points
% dx = 1; 
% dy = 1; 
% 
% scatter3(Z(:,1),Z(:,2),Z(:,3));
% text(Z(:,1)+dx, Z(:,2)+dy, Z(:,3), labels, 'FontSize', 18);

% visualization, courtesy of James Eddy
gene_clustergram = ...
    clustergram(gene_scores, 'ColumnLabels', labels);
permRows = str2num(cell2mat(get(gene_clustergram,'RowLabels')));

close all

% Generate hierarchical clustering distance & linkage to make dendrogram
Y = pdist(gene_scores');
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
imagesc(gene_scores(flipud(permRows),permCols))

% Adjust appearance and add column labels
set(gca,'ytick',[],'xtick',[])
ylabel('Genes','fontsize',8,'fontweight','bold')
for i = 1:12
    text(i,size(gene_scores,1)+48,labels(permCols(i)),...
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