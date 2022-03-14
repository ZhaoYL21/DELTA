
addpath('.\General');
load('.\data\Recon3D_301\Recon3DModel_301.mat');
load('.\data\Recon3D_301\equation.mat');
Recon3DModel_inf = addReaction( Recon3DModel,'biomass_virus',cell2mat(equation));

%%
%Target prediction based on differentially expressed genes obtained from differential expression analysis of transcriptome data
covid_path = 'RESULT_SAVE_PATH_TRANS';
dirs = dir(covid_path);
%读取所有样本目录
file_list={};
for i=1:size(dirs,1)
    isdir=dirs(i).isdir;
    if isdir==1
        name=dirs(i).name;
        if ~isempty(strfind(name,'GSE'))             
            file_list=[file_list;name];
        end
    end
end

if size(file_list,1)==1
    file_list=file_list';
end

for i = 1:size(file_list,1)
    series_name = file_list{i};
    cur_path = strcat(covid_path,series_name,'\');
    disp(cur_path);
    Main_Process( Recon3DModel_inf , cur_path , 'all_Gene-Entrez_up.txt' , 'biomass_virus' , 0 , 1 , 0 , 1 );
end

%%
%Target prediction based on significantly up-regulated genes obtained from the interaction of SARS-CoV-2 proteins with human genes
clear covid_path
covid_path = 'RESULT_SAVE_PATH_INTER';
Main_Process( Recon3DModel_inf , covid_path , 'core_gene_union.txt' , 'biomass_virus' , 0 , 1 , 0 , 1 )

%%
rmpath('.\General');