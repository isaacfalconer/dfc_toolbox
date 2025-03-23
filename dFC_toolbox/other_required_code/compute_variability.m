function dFC_project = compute_variability(dFC_project)

method_choices = {
    'Zhang et al. Formula'
    'Standard deviation'
    };


if ~isfield(dFC_project,'var_method')
    [method,tf] = listdlg('SelectionMode','single','ListString',method_choices);
    if tf==0
        method=1;
    end
else
    method = dFC_project.var_method;
end

% % Need to finish writing this code to exclude lesioned nodes
% cd(CONN_preprocessing_dir)
% cd ..
% pct_spared = readmatrix("pct_spared.csv");

dFC = dFC_project.dFC;
dims = size(dFC);
nroi = dims(1);
nvol = dims(3);

if method==1
    dFC_project.tv_def = 'Zhang et al. Formula';

    corr = zeros(nvol,nvol);
    ncon = numel(dFC_project.cond_names);

    if numel(dims)==3
        var = zeros(nroi,1);
        for j=1:nroi
            for n=1:nvol
                for m=1:nvol
                    if m==n
                        corr(m,n) = NaN;
                    elseif m>n
                        corr_temp = corrcoef(dFC(j,:,m),dFC(j,:,n));
                        corr(m,n) = corr_temp(2,1);
    %                     corr(n,m) = corr_temp(2,1);
                    end
                end
            end
    %         corr_lin = reshape(corr,m*n,1);
    %         var(j) = 1-mean(corr_lin(~isnan(corr_lin)));
            var(j) = 1-mean(corr(~isnan(corr)));
    %         var(isnan(var)) = 0;
            if rem(j,10)==0
                disp([num2str(j),' ROIs done'])
            end
        end
        vars{1} = var;
        mean_var = mean(var,1);
    elseif numel(dims)==4 && ncon==1
        nsub = dims(4);
        var = zeros(nroi,nsub);
        vars = cell(nsub,1);
        mean_var = zeros(nsub,1);
        for i=1:nsub
            disp(['Computing temporal variability for Subject ',num2str(i),' ...'])
            for j=1:nroi
                for n=1:nvol
                    for m=1:nvol
                        if m==n
                            corr(m,n) = NaN;
                        elseif m>n
                            corr_temp = corrcoef(dFC(j,:,m,i),dFC(j,:,n,i));
                            corr(m,n) = corr_temp(2,1);
    %                         corr(n,m) = corr_temp(2,1);
                        end
                    end
                end
    %             corr_lin = reshape(corr,m*n,1);
    %             var(j,i) = 1-mean(corr_lin(~isnan(corr_lin)));
                var(j) = 1-mean(corr(~isnan(corr)));
    %             var(isnan(var)) = 0;
                if rem(j,10)==0
                    disp([num2str(j),' ROIs done'])
                end
            end
            vars{i} = var;
            mean_var(i) = mean(var,1);
        end
    elseif numel(dims)==4 && ncon>1
        var = zeros(nroi,1);
        vars = cell(ncon,1);
        mean_var = zeros(ncon,1);
        for k=1:ncon
            disp(['Computing temporal variability for ',dFC_project.cond_names{k},' ...'])
            for j=1:nroi
                for n=1:nvol
                    for m=1:nvol
                        if m==n
                            corr(m,n) = NaN;
                        elseif m>n
                            corr_temp = corrcoef(dFC(j,:,m,k),dFC(j,:,n,k));
                            corr(m,n) = corr_temp(2,1);
                        end
                    end
                end
    %             corr_lin = reshape(corr,m*n,1);
    %             var(j) = 1-mean(corr_lin(~isnan(corr_lin)));
                var(j) = 1-mean(corr(~isnan(corr)));
    %             var(isnan(var)) = 0;
                if rem(j,10)==0
                    disp([num2str(j),' ROIs done'])
                end
            end
            vars{k} = var;
            mean_var(k) = mean(var,1);
        end
    elseif numel(dims)==5
        nsub = dims(5);
        var = zeros(nroi,1);
        vars = cell(nsub,ncon);
        mean_var = zeros(nsub,ncon);
        for k=1:ncon
            for i=1:nsub
                disp(['Computing temporal variability for Subject ',num2str(i),' ',dFC_project.cond_names{k},' ...'])
                for j=1:nroi
                    for n=1:nvol
                        for m=1:nvol
                            if m==n
                                corr(m,n) = NaN;
                            else
                                corr_temp = corrcoef(dFC(j,:,m,k,i),dFC(j,:,n,k,i));
                                corr(m,n) = corr_temp(2,1);
                            end
                        end
                    end
    %                 corr_lin = reshape(corr,m*n,1);
    %                 var(j) = 1-mean(corr_lin(~isnan(corr_lin)));
                    var(j) = 1-mean(corr(~isnan(corr)));
    %                 var(isnan(var)) = 0;
                    if rem(j,10)==0
                        disp([num2str(j),' ROIs done'])
                    end
                end
                vars{i,k} = var;
                disp(mean(var(~isnan(var))))
                mean_var(i,k) = mean(var(~isnan(var)));
            end
        end
    else
        disp oops!
        return
    end

    dFC_project.vars = vars;
    dFC_project.mean_var = mean_var;
    
elseif method==2
    ncon = numel(dFC_project.cond_names);
    dFC_project.tv_def = 'Std. Dev.';

    if numel(dims)==3
        disp case1
        var = zeros(nroi,nroi);
        roi_var = zeros(nroi,nroi);
        mean_roi_var = zeros(nroi,1);
        for n=1:nroi
            for m=1:nroi
                if m==n
                    var(m,n) = NaN;
                elseif m>n
                    var(m,n) = std(dFC(m,n,:));
                    var(n,m) = var(m,n);
                    roi_var(m,n) = var(m,n);
                end
            end
            roi_var_n = roi_var(:,n);
            mean_roi_var(n) = mean(roi_var_n(~isnan(roi_var_n)));
        end
        vars{1} = mean_roi_var;
        mean_var = mean(var(~isnan(var)));
    elseif numel(dims)==4 && ncon==1
        disp case2
        nsub = dims(4);
        vars = cell(nsub,1);
        mean_var = zeros(nsub,1);
        roi_var = zeros(nroi,nroi);
        mean_roi_var = zeros(nroi,1);
        for i=1:nsub
            var = zeros(nroi,nroi);
            for n=1:nroi
                for m=1:nroi
                    if m==n
                        var(m,n) = NaN;
                    elseif m>n
                        var(m,n) = std(dFC(m,n,:,i));
                        var(n,m) = var(m,n);
                        roi_var(m,n) = var(m,n);
                    end
                end
                roi_var_n = roi_var(:,n);
                mean_roi_var(n) = mean(roi_var_n(~isnan(roi_var_n)));
            end
            vars{i} = mean_roi_var;
            mean_var(i) = mean(var(~isnan(var)));
        end
    elseif numel(dims)==4 && ncon>1
        disp case3
        vars = cell(1,ncon);
        mean_var = zeros(1,ncon);
        roi_var = zeros(nroi,nroi,ncon);
        mean_roi_var = zeros(nroi,ncon);
        for k=1:ncon
            var = zeros(nroi,nroi);
            for n=1:nroi
                for m=1:nroi
                    if m==n
                        var(m,n) = NaN;
                    elseif m>n
                        var(m,n) = std(dFC(m,n,:,k));
                        var(n,m) = var(m,n);
                        roi_var(m,n,k) = var(m,n);
                    end
                end
                roi_var_n = roi_var(~isnan(roi_var(:,n,k)),n,k);
                mean_roi_var(n,k) = mean(roi_var_n,1);
            end
            mean_var(1,k) = mean(var(~isnan(var)));
        end
        vars{1} = mean_roi_var;
    elseif numel(dims)==5
        nsub = dims(5);
        if isfield(dFC_project,'excluded_ROIs') && max(dFC_project.excluded_ROIs>0,[],'all')
            for i=1:nsub
                if dFC_project.which_hemi==1
                    excluded_ROIs = find(dFC_project.excluded_ROIs(:,i));
                    excluded_ROIs = ismember( ...
                        dFC_project.include_hemi_sorted,excluded_ROIs);
                else
                    excluded_ROIs = [];
                    waitfor(warndlg('Exclusion of lesioned ROIs currently not supported for single-hemisphere analyses'))
                end
                dFC(excluded_ROIs,:,:,:,i) = nan;
                dFC(:,excluded_ROIs,:,:,i) = nan;
            end
        end
        dFC_project.dFC_lesion_excluded = dFC;
        numnan = sum(isnan(dFC(:,:,23,1,2)),'all');
        disp(['BU02 number of excluded ROIs: ',num2str(nroi-sqrt((nroi^2)-numnan))]);
        disp(['BU02 expected number of excluded ROIs: ',num2str(sum(ismember( ...
            dFC_project.include_hemi_sorted,find(dFC_project.excluded_ROIs(:,2)))))])
        vars = cell(nsub,1); % cell array to be populated with 
        mean_var = zeros(nsub,ncon);
        roi_var = zeros(nroi,nroi);
        mean_roi_var = zeros(nroi,ncon);
        for i=1:nsub
            for k=1:ncon
                disp(['Computing variability for subject ',num2str(i),' condition ',num2str(k)])
                var = zeros(nroi,nroi);
                for n=1:nroi
                    for m=1:nroi
                        if m==n
                            var(m,n) = NaN;
                        elseif m>n
                            if any(~isnan(dFC(m,n,:,k,i)))
                                var(m,n) = std(dFC(m,n,:,k,i));
                            else
                                var(m,n) = nan;
                            end
                            var(n,m) = var(m,n);
                            roi_var(m,n,k) = var(m,n);
                            roi_var(n,m,k) = var(n,m);
                        end
                    end
                    roi_var_n = roi_var(~isnan(roi_var(:,n,k)),n,k);
                    mean_roi_var(n,k) = mean(roi_var_n,1);
                end
                mean_var(i,k) = mean(var(~isnan(var)));
            end
            vars{i} = mean_roi_var;
            % if i==1
            %     disp(var)
            % end
        end
    else
        disp oops!
        return
    end



    dFC_project.vars = vars;
    dFC_project.mean_var = mean_var;
else
    disp oops!
end

end