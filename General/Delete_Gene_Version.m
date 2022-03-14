function [ Model ] = Delete_Gene_Version( Model )
%DELETE_GENE_VERSION Delete the gene version number of the model

genes = Model.genes;
for i=1:numel(genes)
    x = strsplit(table2array(genes(i,1)),'.');
    genes(i,1) = cellstr(x(1,1));
end
Model.genes = genes;

end

