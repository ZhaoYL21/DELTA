function[ result ] = trans2cell(doublearray)
%Convert double array to cell tuple in str format

result = cell(size(doublearray,1),1);
for i = 1:size(doublearray,1)
    result{i} = num2str(doublearray(i));
end
end