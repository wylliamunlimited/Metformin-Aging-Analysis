% BE 700 A1 Fall 2024
% Final Project, Data Exploration
% Cal Parise, 11/30/2024

brain_data = readtable("Human Brain Data Clean.csv");

gene_id = brain_data.IDENTIFIER;

brain_data_clean = table2array(removevars(brain_data,[1 2]))';
[patients,genes] = size(brain_data_clean);

young = log2(brain_data_clean(1:12,:)); % first 12 patients are young adults (<40yo)
aged = log2(brain_data_clean(13:end,:)); % next are aged as defined in paper


[hvals,pvals] = ttest2(young,aged);
zscores = (-1 * log10(pvals))';
pvals = pvals';

genes_sorted = table(gene_id, pvals, zscores);
genes_sorted = sortrows(genes_sorted,2); % array of gene (probe) names sorted by p-value and z-score

brain_significant_pvals = genes_sorted.pvals(genes_sorted.pvals<=0.05);
brain_significant_genes = genes_sorted.gene_id(genes_sorted.pvals<=0.05);

writecell([brain_significant_genes num2cell(brain_significant_pvals)],"diff_exp_brain_young_vs_aged");

brain_significant_pvals_top1000 = genes_sorted.pvals(1:1000);
brain_significant_genes_top1000 = genes_sorted.gene_id(1:1000);

writecell([brain_significant_genes_top1000 num2cell(brain_significant_pvals_top1000)],"diff_exp_brain_young_vs_aged_top1000");

categories = [];
for i = 1:41
    if i<=12
        categories = [categories,"young"];
    else
        categories = [categories,"aged"];
    end
end

brain_d_z = [brain_data_clean' zscores];
brain_d_z = flipud(sortrows(brain_d_z,patients+1));

new_brain = brain_d_z(1:1000,1:patients)';

[coeff,~,~,~,explained] = pca(brain_data_clean,"Centered",true);

brain_projected = brain_data_clean * coeff;

young_proj = brain_projected(1:12,:);
aged_proj = brain_projected(13:end,:);

figure(2)
hold on
scatter(young_proj(:,1),young_proj(:,2),"bo");
scatter(aged_proj(:,1),aged_proj(:,2),"rs");

title("Principal Component Analysis of Young and Aged Samples","FontWeight","bold")
xlabel("PC1, "+num2str(explained(1))+"%","FontWeight","bold")
ylabel("PC2, "+num2str(explained(2))+"%","FontWeight","bold")
legend("Young","Aged");
grid on
hold off

rng("default");
tsne_brain_new = tsne(new_brain,"Distance","correlation","Standardize",true,"NumDimensions",2,"NumPCAComponents",50);
figure(4);
gscatter(tsne_brain_new(:,1),tsne_brain_new(:,2),categories');
title("t-SNE of Patients, top 1000 differentially-expressed genes","FontSize","bold");