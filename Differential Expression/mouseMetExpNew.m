% BE 700 A1 Fall 2024
% Final Project, Data Exploration
% Cal Parise, 12/12/2024

controlmet_data_table = readtable("mouse_metformin_control_vs_met_expression.tsv","FileType","text","Delimiter","\t");

control_vs_met_table = removevars(controlmet_data_table,["P_Value","t","B","logFC","GI","Gene_title"]);
head(control_vs_met_table)
control_vs_met_p_vals = control_vs_met_table.adj_P_Val(control_vs_met_table.adj_P_Val<=0.05);
control_vs_met_genes = control_vs_met_table.Gene_symbol(control_vs_met_table.adj_P_Val<=0.05);
control_vs_met_ID = control_vs_met_table.ID(control_vs_met_table.adj_P_Val<=0.05);

writecell([control_vs_met_ID control_vs_met_genes num2cell(control_vs_met_p_vals)],"diff_exp_control_vs_met.txt");

controlcr_data_table = readtable("mouse_metformin_control_vs_cr_expression.tsv","FileType","text","Delimiter","\t");

control_vs_cr_table = removevars(controlcr_data_table,["P_Value","t","B","logFC","GI","Gene_title"]);
head(control_vs_cr_table)
control_vs_cr_p_vals = control_vs_cr_table.adj_P_Val(control_vs_cr_table.adj_P_Val<=0.05);
control_vs_cr_genes = control_vs_cr_table.Gene_symbol(control_vs_cr_table.adj_P_Val<=0.05);
control_vs_cr_ID = control_vs_cr_table.ID(control_vs_cr_table.adj_P_Val<=0.05);

writecell([control_vs_cr_ID control_vs_cr_genes num2cell(control_vs_cr_p_vals)],"diff_exp_control_vs_cr.txt");

crmet_data_table = readtable("mouse_metformin_cr_vs_met_expression.tsv","FileType","text","Delimiter","\t");

cr_vs_met_table = removevars(crmet_data_table,["P_Value","t","B","logFC","GI","Gene_title"]);
head(cr_vs_met_table)
cr_vs_met_p_vals = cr_vs_met_table.adj_P_Val(cr_vs_met_table.adj_P_Val<=0.05);
cr_vs_met_genes = cr_vs_met_table.Gene_symbol(cr_vs_met_table.adj_P_Val<=0.05);
cr_vs_met_ID = cr_vs_met_table.ID(cr_vs_met_table.adj_P_Val<=0.05);

writecell([cr_vs_met_ID cr_vs_met_genes num2cell(cr_vs_met_p_vals)],"diff_exp_cr_vs_met.txt");