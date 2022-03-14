function [ ess_gene_name ] = Ess_Gene_Test( Specific_Model , Path_Used , Save_Path , Output_File_Name)
%ESS_GENE_TEST Detection of essential genes

if nargin < 4 || ~exist('Output_File_Name','var')
    Output_File_Name = 'ess_gene_name';
end
if ~exist(Save_Path,'dir')
    mkdir(Save_Path);
end

[ess_gene_name,result,Remove_rxns_result,Remove_rxns_name]=Remove_Gene(Specific_Model , Path_Used);
if size(ess_gene_name,1)==1
    ess_gene_name=ess_gene_name';
end
save(strcat(Save_Path,'\',Output_File_Name,'.mat'),'ess_gene_name');
Save_to_Txt(strcat(Save_Path,'\',Output_File_Name,'.txt'),ess_gene_name);
end

