function dFC_project = transition_freq(dFC_project)

transition_threshold = dFC_project.transition_threshold;

comm_sim = dFC_project.comm_sim;
comm_sim_states = zeros(size(comm_sim));
nsub = size(comm_sim,3);
nwin = size(comm_sim,1);
TF = zeros(nsub,1);

fig = uifigure;
progress_bar = uiprogressdlg(fig,'Title','Computing transition frequency...','Message','');
pause(0.01)
% disp('Computing transition frequency...')

min_state_length = 1; % ceil(5/dFC_project.TR);
% q97 = quantile(comm_sim,0.97,'all');
% dFC_project.transition_threshold = zeros(nsub,1);
% dFC_project = get_transition_threshold(dFC_project);

for j=1:nsub
    % comm_sim_filt(:,:,j) = imgaussfilt(comm_sim(:,:,j),0.7);
    comm_sim_filt(:,:,j) = comm_sim(:,:,j);
    % q95 = quantile(comm_sim,0.95,'all');
    % dFC_project.transition_threshold(j) = q95;

    progress_bar.Message = ['Computing transition frequency for subject ',num2str(j),'...'];
    progress_bar.Value = (j-1)/nsub;
    pause(0.01)

    dFC_sub = dFC_project.dFC(:,:,:,1,j);

    transition = zeros(nwin,1);
    
    transitions = find_transitions(comm_sim_filt(:,:,j), min_state_length, transition_threshold(j));
    % new = 1;
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

    % comm_sim_states(:,:,j) = imgaussfilt(comm_sim(:,:,j),2);
    comm_sim_states(:,:,j) = comm_sim(:,:,j);
    comm_sim_states(:,:,j) = comm_sim_states(:,:,j) - comm_sim_states(:,:,j).*eye(nwin) + eye(nwin);
    comm_sim_states(:,:,j) = comm_sim(:,:,j) - comm_sim(:,:,j).*eye(nwin) + eye(nwin);
    % comm_sim_states(transitions(:,1),:,j) = 1;
    % comm_sim_states(:,transitions(:,1),j) = 1;

    % for i=1:size(transitions,1)
    %         comm_sim_states(transitions(i,1),:,j) = 1;
    %         comm_sim_states(:,transitions(i,1),j) = 1;
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
dFC_project.comm_sim = comm_sim;
dFC_project.comm_sim_states = comm_sim_states;
dFC_project.TF_comm_struct = TF;

% heatmap(w2wcsc_states(:,:,sub)); grid off
% figure(); heatmap(w2wcsc(:,:,sub)); grid off
