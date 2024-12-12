% BE 700 A1 Fall 2024
% Final Project, Data Exploration
% Cal Parise, 11/29/2024

clear vars, close all

diab_geneexp = readtable("Insulin Resistance Gene Expression.csv");
diab_pheno = readtable("Insulin Resistance Phenotypes.csv");

gene_labels = diab_geneexp.Var1;
si = diab_pheno.si_ch1;
diab_geneexp_clean = table2array(removevars(diab_geneexp,"Var1"));
[genes,patients] = size(diab_geneexp_clean);

ID_a1c = table(diab_pheno.Var1, diab_pheno.hemoglobinA1c_ch1);
head(ID_a1c)

diab_categories = zeros(patients,1);
for i = 1:patients

    if ID_a1c.Var2(i) < 5.7
        diab_categories(i) = 0;
    else
        diab_categories(i) = 1;
    end

end

ID_a1c.Category = diab_categories;
head(ID_a1c)

diab_geneexp_clean = [diab_geneexp_clean' ID_a1c.Category];

geneexp_sorted = sortrows(diab_geneexp_clean,genes+1);

normal_a1c = geneexp_sorted(geneexp_sorted(:,end)==0,1:genes);
prediabetic_and_diabetic = geneexp_sorted(geneexp_sorted(:,end)==1,1:genes);

[hvals,pvals] = ttest2(normal_a1c,prediabetic_and_diabetic);
zscores = (-1 * log10(pvals))';
pvals = pvals';

genes_sorted = table(gene_labels, pvals, zscores);
genes_sorted = sortrows(genes_sorted,2); % array of gene (probe) names sorted by p-value and z-score

human_muscle_significant_pvals = genes_sorted.pvals(genes_sorted.pvals<=0.05);
human_muscle_significant_genes = genes_sorted.gene_labels(genes_sorted.pvals<=0.05);

writecell([human_muscle_significant_genes num2cell(human_muscle_significant_pvals)],"diff_exp_human_muscle_normal_vs_diabetic_a1c");