function [ ] = Main_Process( Background_Model , Path , Input_File_Name , Target_Object , Need_Map_Gup , Need_Map_Toxicity , Need_Split_Gup , Need_Split_Toxicity )
%MAIN_PROCESS 

if nargin < 8 || ~exist('Need_Split_Toxicity','var')
    Need_Split_Toxicity = 1;
end
if nargin < 7 || ~exist('Need_Split_Gup','var')
    Need_Split_Gup = 0;
end
if nargin < 6 || ~exist('Need_Map_Toxicity','var')
    Need_Map_Toxicity = 1;
end
if nargin < 5 || ~exist('Need_Map_Gup','var')
    Need_Map_Gup = 0;
end

Path_Used = Create_Path_Used();

Background_Model_splite = Delete_Gene_Version(Background_Model);

Map=struct();
Map.bool=0;
if Need_Map_Gup == 1
   [Entrez,ENSG] = xlsread(Path_Used.map_path); 
   Map.bool=1;
   Map.Entrez=Entrez;
   Map.ENSG=ENSG;
end
Gene_up=Select_Gene_Up(Background_Model_splite , Path , Input_File_Name , Path , 'Gene_up' , Map , Need_Split_Gup);

Reac_up=Map_Gene_to_Rxn( Gene_up , Background_Model_splite , Target_Object , Path , 'Reac_up' );

Specific_Model=Reconstruction( Reac_up , Background_Model , Target_Object , Path_Used , Path );

ess_gene_name=Ess_Gene_Test( Specific_Model , Path_Used , Path , 'ess_gene_name');

Map=struct();
Map.bool=0;
if Need_Map_Toxicity == 1
   [Entrez,ENSG] = xlsread(Path_Used.map_path); 
   Map.bool=1;
   Map.Entrez=Entrez;
   Map.ENSG=ENSG;
end
target=Toxicity_Test( ess_gene_name , Path_Used , Path , 'Toxicity' , 'Target' , Map , Need_Split_Toxicity );

end

