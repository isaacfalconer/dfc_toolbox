function transitions = find_transitions(A, min_state_length, transition_threshold)

    nwin = size(A,1);
    if nwin < 2*min_state_length
        errordlg('Error using ''transitions'' function: duration must be at least 2 times the minimum state length')
    end

    transitions = zeros(0,1);
    last = nwin-2*min_state_length-1;
    % for i=2:last
        % stop = i+min_state_length-1;
        % C = A((i+min_state_length):(stop+min_state_length),i:stop);
        % a = mean(C(C>0));
        % if a<transition_threshold
            % transitions(size(transitions,1)+1,1:2) = [stop a];
    for i=3:nwin
        if A(i,i-1)<transition_threshold
            transitions(size(transitions,1)+1,1) = i;
        elseif A(i,i-2)<transition_threshold
            transitions(size(transitions,1)+1,1) = i-1;
        end
    end

    if size(transitions,1)==0
        return
    end
    
    transitions = unique(transitions);

    transitions = transitions(transitions(:,1)>0,:);
    % transitions(:,3) = 0;
    % transitions = transitions(transitions(:,2)==min(transitions(:,2)),:);
    % if size(transitions,1) > 1
    %     transitions = transitions(1,:);
    % end
end