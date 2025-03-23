function [comm_diff, perm_needed, V1, V2] = find_comm_diff(V1, V2)

comm_diff1 = sum(V1~=V2)/length(V1);

count = 0;
for i=1:max(V1)
    if sum(V1==i)<2
        V1(V1==i) = 0;
        count = count + 1;
    else
        V1(V1==i) = i-count;
    end
end

count = 0;
for i=1:max(V2)
    if sum(V2==i)<2
        V2(V2==i) = 0;
        count = count + 1;
    else
        V2(V2==i) = i-count;
    end
end

if max(V2)>8
    p = zeros(1000,max(V2));
    for i=1:4000
        next = randperm(max(V2));
        p(i,:) = next;
    end
else
    p = perms(1:max(V2));
end

diffs = zeros(max(V2),1);

for i=1:size(p,1)
    V3 = zeros(size(V2));
    for j=1:max(V2)
        V3(V2==j) = p(i,j);
    end
    diffs(i) = sum(V1~=V3)/length(V1);
end

max_diff = (length(V1)/max(V1)-1)/(length(V1)/max(V1));

comm_diff = min(diffs); % /max_diff;

perm_needed = comm_diff < comm_diff1;

if ~perm_needed
    comm_diff = comm_diff1;
end

end