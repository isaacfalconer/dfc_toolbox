function [TF,transition] = compute_TF(sub_clust,TR)

nsub = size(sub_clust,2);
nvol = size(sub_clust,1);
duration = nvol*TR/60;

TF = zeros(nsub,1);

transition = zeros(nvol,nsub);
for i=1:nsub
    transitions = zeros(nvol,1);
    for j=2:nvol
        if sub_clust(j,i) ~= sub_clust(j-1,i)
            transitions(j) = 1;
        end
    end
    TF(i) = sum(transitions)/duration;
    transition(:,i) = transitions;
end

end