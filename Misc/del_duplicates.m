function M = del_duplicates(data)
    %Delete duplicates rows in data matrix while preserving original order
    if size(data,1)>1 && size(data,2)>1
        [M,ind] = unique(data, 'rows', 'first');
        [~,ind] = sort(ind);
        M = M(ind,:);
    else
        [M,ind] = unique(data, 'first');
        [~,ind] = sort(ind);
        M = M(ind);
    end
end