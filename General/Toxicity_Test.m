function [ target ] = Toxicity_Test( ess_gene_name , Path_Used , Save_Path , Output_File_Name , Target_File_Name , Need_Map , Need_Split)
%TOXICITY_TEST Toxicity testing to predict candidate targets

if nargin < 7 || ~exist('Need_Split','var')
    Need_Split = 1;
end
if nargin < 6 || ~exist('Need_Map','var')
    Need_Map = struct();
    Need_Map.bool = 0;
end
if nargin < 5 || ~exist('Target_File_Name','var')
    Target_File_Name = 'Target';
end
if nargin < 4 || ~exist('Output_File_Name','var')
    Output_File_Name = 'Toxicity';
end
if ~exist(Save_Path,'dir')
    mkdir(Save_Path);
end

addpath(Path_Used.healthy_path);

all_ess_gene=ess_gene_name;
if size(all_ess_gene,1)==1
    all_ess_gene=all_ess_gene';
end
gene_size=size(all_ess_gene,1);
addpath(Path_Used.raven_core_path);
addpath(Path_Used.raven_solver_path);
task=parseTaskList('TASKS.xlsx');%56 metabolic tasks
%83 normal models
datafile={'vagina - squamous epithelial cells', 'appendix - glandular cells', 'appendix - lymphoid tissue', 'bone marrow - hematopoietic cells', 'breast - adipocytes', 'breast - glandular cells','breast - myoepithelial cells', 'bronchus - respiratory epithelial cells', 'cerebellum - cells in granular layer', 'cerebellum - cells in molecular layer', 'cerebellum - Purkinje cells', 'cerebral cortex - endothelial cells','cerebral cortex - glial cells', 'cerebral cortex - neuronal cells', 'cerebral cortex - neuropil', 'cervix, uterine - glandular cells', 'cervix, uterine - squamous epithelial cells', 'colon - endothelial cells','colon - glandular cells', 'colon - peripheral nerve-ganglion', 'duodenum - glandular cells', 'epididymis - glandular cells', 'esophagus - squamous epithelial cells', 'fallopian tube - glandular cells','gallbladder - glandular cells', 'heart muscle - myocytes', 'hippocampus - glial cells', 'hippocampus - neuronal cells', 'kidney - cells in glomeruli', 'kidney - cells in tubules','lateral ventricle - glial cells', 'lateral ventricle - neuronal cells', 'liver - bile duct cells', 'liver - hepatocytes', 'lung - macrophages', 'lung - pneumocytes','lymph node - germinal center cells', 'lymph node - non-germinal center cells', 'nasopharynx - respiratory epithelial cells', 'oral mucosa - squamous epithelial cells', 'ovary - follicle cells', 'ovary - ovarian stroma cells','pancreas - exocrine glandular cells', 'pancreas - islets of Langerhans', 'parathyroid gland - glandular cells', 'placenta - decidual cells', 'placenta - trophoblastic cells', 'prostate - glandular cells','rectum - glandular cells', 'salivary gland - glandular cells', 'seminal vesicle - glandular cells', 'skeletal muscle - myocytes', 'skin 1 - fibroblasts', 'skin 1 - keratinocytes','skin 1 - Langerhans', 'skin 1 - melanocytes', 'skin 2 - epidermal cells', 'small intestine - glandular cells', 'smooth muscle - smooth muscle cells', 'soft tissue 1 - adipocytes','soft tissue 1 - chondrocytes', 'soft tissue 1 - fibroblasts', 'soft tissue 1 - peripheral nerve', 'soft tissue 2 - adipocytes', 'soft tissue 2 - chondrocytes', 'soft tissue 2 - fibroblasts','soft tissue 2 - peripheral nerve', 'spleen - cells in red pulp', 'spleen - cells in white pulp', 'stomach 1 - glandular cells', 'stomach 2 - glandular cells', 'testis - cells in seminiferous ducts','testis - Leydig cells', 'thyroid gland - glandular cells', 'tonsil - germinal center cells', 'tonsil - non-germinal center cells', 'tonsil - squamous epithelial cells', 'urinary bladder - urothelial cells','uterus 1 - cells in endometrial stroma', 'uterus 1 - glandular cells', 'uterus 2 - cells in endometrial stroma', 'uterus 2 - glandular cells', 'vagina - squamous epithelial cells'};
type='.mat';
m_file=strcat(datafile,type);
fail_task=cell(gene_size,83);
result_gene=83*ones(size(all_ess_gene,1),1);
flag_gene=zeros(size(all_ess_gene,1),83);

rmpath(Path_Used.cobra_refine_path);
if Need_Map.bool==1
    ENSG=Need_Map.ENSG;
    Entrez=Need_Map.Entrez;
end
load(Path_Used.acknow_path);

for j=1:83
    
    load(datafile{j});
    for i=1:size(all_ess_gene,1)
        Rgene=all_ess_gene{i};
        if Need_Split==1
            Rgene_split=strsplit(Rgene,'.');
            Rgene=Rgene_split{1};
        end
        if Need_Map.bool==1
            Rgene=str2num(Rgene);
            [~,loc]=ismember(Rgene,Entrez);
            data_loc=find(loc>0);
            loc=loc(data_loc);
            Rgene=ENSG(loc);
        end
        [Rtf,Rloc]=ismember(Rgene,Acknow_Gene);
        if(Rtf==0)
            if(ismember(Rgene,model.genes))
                Model_reduced=removeGenes(model,Rgene);
                taskReport=checkTasks(Model_reduced,[],true,false,false,task);
                t=find(taskReport.ok==0);
                if t
                    result_gene(i,1)=result_gene(i,1)-1;
                    flag_gene(i,j)=-1;
                    fail_task(i,j)={t};
                else
                    flag_gene(i,j)=1;
                end
            else
                result_gene(i,1)=result_gene(i,1)-1;
            end
            if(j==83)
                Acknow_Gene{size(Acknow_Gene,2)+1}=Rgene{1};
                Acknow_Num=[Acknow_Num;result_gene(i,1)];
            end
        else
            result_gene(i,1)=Acknow_Num(Rloc);
        end
    end  
    
end

save(Path_Used.acknow_path,'Acknow_Gene','Acknow_Num');
save(strcat(Save_Path,'\',Output_File_Name,'.mat'),'result_gene','fail_task','flag_gene');
target=all_ess_gene(result_gene>42);
save(strcat(Save_Path,'\',Target_File_Name,'.mat'),'target');
Save_to_Txt(strcat(Save_Path,'\',Target_File_Name,'.txt'),target);

rmpath(Path_Used.healthy_path);
rmpath(Path_Used.raven_core_path);
rmpath(Path_Used.raven_solver_path);
addpath(Path_Used.cobra_refine_path);

end

