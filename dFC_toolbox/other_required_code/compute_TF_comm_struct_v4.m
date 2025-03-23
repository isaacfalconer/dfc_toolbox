function dFC_project = compute_TF_comm_struct_v4(dFC_project, condition)


% condition = 1;

nsub = size(dFC_project.dFC,5);
% nsub = 1;

TF = zeros(nsub,1);
mean_duration = zeros(nsub,1);
time_in_state = zeros(nsub,1);
mean_in_state_dFC = zeros(nsub,1);
mean_out_of_state_dFC = zeros(nsub,1);
mean_dFC = zeros(nsub,1);

nroi = size(dFC_project.dFC,1);
nwin = size(dFC_project.dFC,3);
comm_sim = zeros(nwin,nwin,nsub); % Window-to-window community structure correlation
comm_sim_filt = zeros(nwin,nwin,nsub);
comm_sim_states = zeros(nwin,nwin,nsub);

dFC_for_kmeans = zeros(nroi,nroi,0);
subs_for_kmeans = zeros(0);
dFC_flat = zeros(0,nroi^2);

trans = cell(nsub,1);
rng_state = rng();

fig = uifigure;
progress_bar = uiprogressdlg(fig,'Title','Computing inter-window similarity...','Message','');
pause(0.01)

dFC_comm_struct = zeros(nroi,nroi,nwin,nsub);
mean_similarity = zeros(nsub,1);
spared = ~any(dFC_project.excluded_ROIs(dFC_project.include_hemi_sorted,:),2);
dFC_project.spared = find(spared);

perm_needed_tot = nan(nsub,1);
dFC_project.spared = cell(nsub,1);
for j=1:nsub
    
    % spared = ~dFC_project.excluded_ROIs(dFC_project.include_hemi_sorted,j);
    % dFC_project.spared{j} = find(spared);

    % disp(['Computing inter-window correlations for subject ',num2str(j),'...'])
    progress_bar.Message = ['Computing inter-window similarity for subject ',num2str(j),'...'];
    progress_bar.Value = (j-1)/nsub;
    pause(0.01)
    dFC_sub = dFC_project.dFC(:,:,:,1,j);
    if dFC_project.which_hemi==1 && any(dFC_project.excluded_ROIs,'all')
        excluded_ROIs = find(dFC_project.excluded_ROIs(:,j));
        excluded_ROIs = ismember( ...
            dFC_project.include_hemi_sorted,excluded_ROIs);
    else
        excluded_ROIs = [];
    end
    dFC_sub(excluded_ROIs,:,:) = 0;
    dFC_sub(:,excluded_ROIs,:) = 0;
    dFC_comm_struct(:,:,:,j) = dFC_sub;
        
    comm = zeros(nroi,nwin);
    for i=1:nwin
        A = dFC_sub(:,:,i);
        % A(A<0) = 0;
        A = A - min(A,[],'all');
        rng(103)
        [comm(:,i),~] = community_louvain(A);
    end

    % if isfield(dFC_project,'w2wcsc')
    %     w2wcsc = dFC_project.w2wcsc;
    %     for i=1:nwin
    %         w2wcsc(i,i,j) = 0;
    %     end
    % else
    perm_needed_tot(j) = 0;
    for n=1:nwin
        for m=n+1:nwin
            comm1 = comm(spared,n);
            comm2 = comm(spared,m);
            % comm1 = comm(:,n);
            % comm2 = comm(:,m);
            [diff, perm_needed] = find_comm_diff(comm1,comm2);
            comm_sim(n,m,j) = 1-diff;
            comm_sim(m,n,j) = 1-diff;
            if perm_needed && m-n<12
                perm_needed_tot(j) = perm_needed_tot(j) + 1;
            end
        end
    end
        comm_sim(1,1,j) = comm_sim(1,2,j);
        comm_sim(nwin,nwin,j) = comm_sim(nwin-1,nwin,j);
        for i=2:nwin-1
            comm_sim(i,i,j) = mean([comm_sim(i-1,i,j),comm_sim(i,i+1,j)]);
        end
        % comm_sim(:,:,j) = imgaussfilt(comm_sim(:,:,j),0.7);
    for i=1:nwin
        comm_sim(i,i,j) = nan;
    end
    comm_sim_j = comm_sim(:,:,j);
    mean_similarity(j) = mean(comm_sim_j(~isnan(comm_sim_j)),'all');
end

dFC_project.comm_sim = comm_sim;

close(progress_bar)
close(fig)

fig = uifigure;
progress_bar = uiprogressdlg(fig,'Title','Computing community stability...','Message','');
pause(0.01)
% disp('Computing transition frequency...')

min_state_length = 1; % ceil(5/dFC_project.TR);
% q97 = quantile(comm_sim,0.97,'all');
dFC_project.transition_threshold = zeros(nsub,1);
% dFC_project = get_transition_threshold(dFC_project);

for j=1:nsub
    % comm_sim_filt(:,:,j) = imgaussfilt(comm_sim(:,:,j),0.7);
    comm_sim_filt(:,:,j) = comm_sim(:,:,j);
    % q95 = quantile(comm_sim,0.95,'all');
    % dFC_project.transition_threshold(j) = q95;

    progress_bar.Message = ['Computing community stability for subject ',num2str(j),'...'];
    progress_bar.Value = (j-1)/nsub;
    pause(0.01)

    dFC_sub = dFC_project.dFC(:,:,:,1,j);

    transition = zeros(nwin,1);
    
    try
        dFC_project.transition_threshold(j) = get_transition_threshold(comm_sim_filt,dFC_project.TR,0);
    catch
        dFC_project.transition_threshold(j) = 1;
    end
    transitions = find_transitions(comm_sim_filt(:,:,j), min_state_length, dFC_project.transition_threshold(j));
    new = 1;
    % while new==1
    %     new_transitions = zeros(0,3);
    %     try
    %         B = comm_sim_filt(1:transitions(1,1),1:transitions(1,1),j);
    %     catch
    %         transitions(1,1)
    %         size(comm_sim_filt(1:transitions(1,1),1:transitions(1,1),j))
    %         return
    %     end
    %     if size(B,1)> 2*min_state_length
    %         new_transitions_first = find_transitions(B, min_state_length, dFC_project.transition_threshold(j));
    %     else
    %         new_transitions_first = zeros(0,4);
    %     end
    %     new_transitions = [new_transitions_first];
    %     last = size(transitions,1);
    %     B = comm_sim_filt(transitions(last,1):nwin,transitions(last,1):nwin,j);
    %     if size(B,1)> 2*min_state_length
    %         new_transitions_last = find_transitions(B, min_state_length, dFC_project.transition_threshold(j)) + [transitions(last,1) 0 0];
    %     else
    %         new_transitions_last = zeros(0,4);
    %     end
    %     new_transitions = [new_transitions; new_transitions_last];
    %     for i=1:size(transitions,1)-1
    %         if transitions(i+1,1) - transitions(i,1) - 1 > 2*min_state_length
    %             B = comm_sim_filt(transitions(i,1)+1:transitions(i+1,1)-1,transitions(i,1)+1:transitions(i+1,1)-1,j);
    %             new_transitions_i = find_transitions(B, min_state_length, dFC_project.transition_threshold(j)) + [transitions(i,1)+1 0 0];
    %             new_transitions = [new_transitions; new_transitions_i];
    %         end
    %     end
    %     new = size(new_transitions,1)>0;
    %     if new
    %         transitions = [transitions; new_transitions];
    %     end
    %     [~,index] = sort(transitions(:,1));
    %     transitions = transitions(index,:);
    % end
            
    transitions = transitions(transitions(:,1)>0,:);

    % for i=1:size(transitions,1)
    %     A(states(i,1)-1:states(i,2),states(i,1)-1) = 2;
    %     A(states(i,1)-1,states(i,1)-1:states(i,2)) = 2;
    %     A(states(i,1)-1:states(i,2),states(i,2)) = 2;
    %     A(states(i,2),states(i,1)-1:states(i,2)) = 2;
    %     if i<size(states,1) && states(i+1,1) - states(i,2) > min_state_length
    %         A(states(i,2)+1:states(i+1,1)-2,states(i,2)+1) = nan;
    %         A(states(i,2)+1,states(i,2)+1:states(i+1,1)-2) = nan;
    %         A(states(i,2)+1:states(i+1,1)-2,states(i+1,1)-2) = nan;
    %         A(states(i+1,1)-2,states(i,2)+1:states(i+1,1)-2) = nan;
    %         for n=states(i,2)+2:states(i+1,1)-3
    %             A(n,n) = comm_sim(n,n,j);
    %         end
    %     elseif i<size(states,1)
    %         for n=states(i,2)+1:states(i+1,1)-1
    %             A(n,n) = nan;
    %         end
    %     end
    % end

    % transition(p<0.05&p2>=0.05) = 1;
    TF(j) = size(transitions,1);

    % for i=1:size(states,1)-1
    %     if states(i+1,1) - states(i,2) > min_state_length
    %         TF(j) = TF(j) + 1;
    %     end
    % end

    TF(j) = TF(j)/(nwin*dFC_project.TR/60);

    % for i=2:nwin
    %     if transition(i)==1
    %         A(i,i) = nan;
    %     end
    % end

    trans{j} = transitions(:,1);
    disp(['Subject ',num2str(j),' done...'])

    % comm_sim_states(:,:,j) = imgaussfilt(comm_sim(:,:,j),0.7);
    comm_sim_states(:,:,j) = comm_sim(:,:,j);
    comm_sim_states(:,:,j) = comm_sim_states(:,:,j) - comm_sim_states(:,:,j).*eye(nwin) + eye(nwin);
    comm_sim_states(:,:,j) = comm_sim(:,:,j) - comm_sim(:,:,j).*eye(nwin) + eye(nwin);
    
    % for i=1:size(transitions,1)
    %     if i==1 || i==size(transitions,1)
    %         comm_sim_states(transitions(i,1),:,j) = 1;
    %         comm_sim_states(:,transitions(i,1),j) = 1;
    %         plot_line = 0;
    %     end
    %     if plot_line
    %         comm_sim_states(transitions(i,1),:,j) = 1;
    %         comm_sim_states(:,transitions(i,1),j) = 1;
    %         plot_line = 0;
    %     end
    %     if i<size(transitions,1)-1 && transitions(i+2) - transitions(i+1) > 1
    %         plot_line = 1;
    %     end
    % 
    % end

    for i=1:size(transitions,1)
        if i==1 || i==size(transitions,1)
            comm_sim_states(transitions(i,1),:,j) = 0;
            comm_sim_states(:,transitions(i,1),j) = 0;
            plot_line = 0;
        else
            if plot_line || transitions(i) - transitions(i-1) > 1 && transitions(i+1) - transitions(i) > 1
                comm_sim_states(transitions(i,1),:,j) = 0;
                comm_sim_states(:,transitions(i,1),j) = 0;
                plot_line = 0;
            end
            if i<size(transitions,1)-1 && transitions(i+2) - transitions(i+1) > 1
                plot_line = 1;
            end
        end
    end

    % for i=1:nwin
    %     for k=1:nwin
    %         % if isnan(A(i,k))
    %         %     comm_sim_states(i,k,j) = 0;
    %         % end
    %         if A(i,k)==2
    %             comm_sim_states(i,k,j) = 1;
    %         end
    %     end
    % end

    % duration = states(:,2) - states(:,1);
    % mean_duration(j) = mean(duration)*dFC_project.TR;

    % dFC_flats = zeros(0,nroi^2);
    % for i=1:size(states,1)
    %     start = states(i,1);
    %     stop = states(i,2);
    %     new_dFC = dFC_sub(:,:,start:stop);
    %     new_dFC_flat = reshape(permute(dFC_sub(:,:,start:stop),[3 1 2]),[stop-start+1 nroi^2]);
    %     dFC_for_kmeans = cat(3, dFC_for_kmeans, new_dFC);
    %     dFC_flats = [dFC_flats; new_dFC_flat];
    % end
    % dFC_flat = [dFC_flat; dFC_flats];
    % subs = j*ones(size(dFC_flats,1),1);
    % subs_for_kmeans = [subs_for_kmeans; subs];

    % comm_sim(:,:,j) = imgaussfilt(comm_sim(:,:,j),2);
    % comm_sim(1,1,j) = nan;
end

close(progress_bar)
close(fig)

% figure; heatmap(comm_sim_states(:,:,1)), grid off

% dFC_project.mean_duration = mean_duration;
% dFC_project.time_in_state = time_in_state;
% dFC_project.mean_in_state_dFC = mean_in_state_dFC;
% dFC_project.mean_out_of_state_dFC = mean_out_of_state_dFC;
% dFC_project.dFC_flat = dFC_flat;
% dFC_project.dFC_for_kmeans = dFC_for_kmeans;
% dFC_project.subs_for_kmeans = subs_for_kmeans;
% dFC_project.perm_needed_tot = perm_needed_tot;
% 
dFC_project.trans = trans;
dFC_project.comm_sim_states = comm_sim_states;
dFC_project.TF_comm_struct = TF;
dFC_project.dFC_comm_struct = dFC_comm_struct;
dFC_project.mean_similarity = mean_similarity;

% heatmap(w2wcsc_states(:,:,sub)); grid off
% figure(); heatmap(w2wcsc(:,:,sub)); grid off

rng(rng_state)

end