function new_tablue = pivot_table(row , col , tablue)

[R,L] = size(tablue);
if tablue(row,col) ~=0
    tablue(row,:) = tablue(row,:)/tablue(row,col);
else
    disp('The Pivot Element can not be zero');
    return
end
for r = 1:1:R
    if r ~= row
        tablue(r,:) = tablue(r,:)-(tablue(r,col)/tablue(row,col))*tablue(row,:);
    end
end
new_tablue = tablue;
end