function [ess_gene_name,result,Remove_rxns_result,Remove_rxns_name]=Remove_Gene(Model , Path_Used)
% Remove_Gene£ºGene knockout is performed to obtain the gene whose objective function is 0 after knockout.

addpath(Path_Used.raven_core_path);
addpath(Path_Used.raven_solver_path);
rmpath(Path_Used.cobra_refine_path);

result=zeros(size(Model.genes,1),1);
Remove_rxns_result=cell(size(Model.genes,1),1);
Remove_rxns_name=cell(size(Model.genes,1),1);

for i=1:size(Model.genes,1)
    Rgene=Model.genes(i,1);
    [Model_reduced,notDeleted,remove_reaction]=removeGenes(Model,Rgene);
    solution=solveLPR(Model_reduced);
    result(i,1)=solution.f;
    Remove_rxns_result{i,1}=remove_reaction;
    Remove_rxns_name{i,1}=Model.rxns(remove_reaction);

end
rmpath(Path_Used.raven_core_path);
rmpath(Path_Used.raven_solver_path);
addpath(Path_Used.cobra_refine_path);
result=-result;
ess_gene1=find(result<1e-6);
ess_gene_name=Model.genes(ess_gene1);
Remove_rxns_result=Remove_rxns_result(ess_gene1);
Remove_rxns_name= Remove_rxns_name(ess_gene1);