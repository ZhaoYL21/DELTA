function [ Specific_Model ] = Reconstruction( Reac_up , Background_Model , Target_Object , Path_Used , Save_Path , Output_File_Name , Output_File_Name_Gene , Output_File_Name_Rxn , Output_File_Name_Met)
%RECONSTRUCTION Generate specific metabolic network models

if nargin < 9 || ~exist('Output_File_Name_Met','var')
    Output_File_Name_Met = 'Model_Met';
end
if nargin < 8 || ~exist('Output_File_Name_Rxn','var')
    Output_File_Name_Rxn = 'Model_Reac';
end
if nargin < 7 || ~exist('Output_File_Name_Gene','var')
    Output_File_Name_Gene = 'Model_Gene';
end
if nargin < 6 || ~exist('Output_File_Name','var')
    Output_File_Name = 'Specific_Model';
end
if ~exist(Save_Path,'dir')
    mkdir(Save_Path);
end


changeCobraSolver('ibm_cplex', 'LP', 0);
addpath(Path_Used.cobra_refine_path);

Specific_Model=fastcore(Background_Model,Reac_up,1e-4,0);
[~,loc]=ismember(Specific_Model.rxns,Target_Object);
Specific_Model.c(find(loc>0))=1;

genes=Specific_Model.genes;
reacs=Specific_Model.rxns;
mets=Specific_Model.metNames;
Save_to_Txt(strcat(Save_Path,'\',Output_File_Name_Gene,'.txt'),genes);
Save_to_Txt(strcat(Save_Path,'\',Output_File_Name_Rxn,'.txt'),reacs);
Save_to_Txt(strcat(Save_Path,'\',Output_File_Name_Met,'.txt'),mets);

save(strcat(Save_Path,'\',Output_File_Name,'.mat'),'Specific_Model');

end

