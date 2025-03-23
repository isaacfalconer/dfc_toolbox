function [threshold,AUC,density_1] = get_transition_threshold(comm_sim,TR,plot_density)

nwin = size(comm_sim,1);
nsub = size(comm_sim,3);

% Keep only similarity values for windows more than 30 seconds apart
comm_sim_1 = zeros(nwin,nwin,nsub);
for i=1:nwin
    for j=1:nwin
        if abs(j-i)<=(60/TR)
            comm_sim_1(i,j,:) = nan;
        else
            comm_sim_1(i,j,:) = comm_sim(i,j);
        end
    end
end

% Keep only similarity values for windows within 30 seconds of each other
comm_sim_2 = zeros(nwin,nwin,nsub);
for i=1:nwin
    for j=1:nwin
        if abs(j-i)>1 || i==j
            comm_sim_2(i,j,:) = nan;
        else
            comm_sim_2(i,j,:) = comm_sim(i,j);
        end
    end
end

pts = (1:100)/100;
[density_1,x1] = ksdensity(reshape(comm_sim_1,[nsub*nwin^2 1]),pts);
[density_2,x2] = ksdensity(reshape(comm_sim_2,[nsub*nwin^2 1]),pts);
density_1 = smoothdata(density_1);
density_2 = smoothdata(density_2);

start = find(density_1 == max(density_1));
stop = find(density_2 == max(density_2));

abs_diff = abs(density_2(start:stop) - density_1(start:stop));
stop = find(abs_diff == min(abs_diff))+start-1;
threshold = pts(stop);

% for i=1:numel(pts)-1
%     AUC = (pts(2)-pts(1))*sum(density_1(numel(pts)-i:numel(pts)));
%     if AUC>=0.05
%         stop = numel(pts)-i;
%         threshold = pts(stop);
%         break
%     end
% end

% nonzero1 = find(density_1>0.001);
% threshold = pts(max(nonzero1));

if plot_density
    figure; hold on
    plot(x1,density_1/100)
    plot(x2,density_2/100)
    ylim([0 1.1*max(density_1/100)])
    xlim([0 1])
    plot([threshold threshold],[0 1.1*max(density_1/100)],'k','LineStyle','--')
    area([1:stop]*(pts(2)-pts(1)),density_2(1:stop)/100,'EdgeColor','none','FaceColor',[0.8500, 0.3250, 0.0980],'FaceAlpha',0.5)
    hold off
    legend({'Windows >60s apart','Adjacent windows','Transition Threshold'},'Location','northeast')
    xlabel('Inter-window similarity')
    ylabel('Probability')
end

% AUC = (pts(2)-pts(1))*sum(density_2(1:intersection));

end
