function [ Gene_up ] = Select_Gene_Up( Background_Model , Load_Path , Input_File_Name , Save_Path , Output_File_Name , Need_Map , Need_Split )
%SELECT_GENE_UP Metabolism-related genes are extracted from a gene set

if nargin < 7 || ~exist('Need_Split','var')
    Need_Split = 0;
end
if nargin < 6 || ~exist('Need_Map','var')
    Need_Map = struct();
    Need_Map.bool = 0;
end
if nargin < 5 || ~exist('Output_File_Name','var')
    Output_File_Name = 'Gene_up';
end
if ~exist(Save_Path,'dir')
    mkdir(Save_Path);
end

up_file = strcat(Load_Path,'\',Input_File_Name);
fid=fopen(up_file);
tline = fgetl(fid);

up_gene = {};

while(ischar(tline))
    str=tline;
    if Need_Split==1
        str=regexp( str, '\w+',  'match' );
        str=str(1);
    end
    up_gene=[up_gene;str];
    tline=fgetl(fid);
end
fclose(fid);
genes=Background_Model.genes;
if Need_Map.bool==1
    ENSG=Need_Map.ENSG;
    Entrez=Need_Map.Entrez;
    [~,loc]=ismember(up_gene,ENSG);
    r=find(loc>0);
    up_gene=Entrez(loc(r));
    up_gene=trans2cell(up_gene);
    clear r
end

[tf,~] = ismember(up_gene,genes);
r=find(tf>0);
Gene_up = up_gene(r);

save(strcat(Save_Path,'\',Output_File_Name,'.mat'),'Gene_up');
Save_to_Txt(strcat(Save_Path,'\',Output_File_Name,'.txt'),Gene_up);
end

