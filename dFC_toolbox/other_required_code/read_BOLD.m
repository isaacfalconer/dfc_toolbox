function [BOLD,ROIs,nvol,dFC_project] = read_BOLD(dFC_project)

%Extracts time series data from .mat files in the directory specified by
%'path' and returns a 4-dimensional matrix of time series data. Excludes 
%confound regressors and GM, WM, and CSF ROIs. Trims time series to the
%length of the shortest scan in the dataset for consitency

%Input Arguments:
%path = directory containing CONN output (filename: 'ROI_SubjectXXX_ConditionXXX.mat')
%subj_num = number of subjects
%cond_num = number of conditions

%Output:
%BOLD = matrix with dimensions [volume number, ROI number, subject number, condition number]

oldPath = pwd;
cd(dFC_project.CONN_preprocessing_dir)

% Get number of conditions and subjects
if ~isfield(dFC_project,'nsub')
    dFC_project = get_ncon_nsub(dFC_project);
end
ncon = dFC_project.ncon;
nsub = dFC_project.nsub;

ncon_actual = numel(dFC_project.conds); % number of conditions to be included in analysis
nvols = zeros(nsub,ncon_actual); % empty matrix to hold the total number of volumes for each subject for each file
nvols_actual = zeros(nsub,ncon_actual); % emtpy matrix to hold number of volumes corresponding to each condition for each subject

% Check if there are files to convert then get time series parameters
if ncon < 1
    disp('No files to convert')
    return
else
    full_roi = numel(load('ROI_Subject001_Condition001.mat').names);
    nroi = dFC_project.nroi; % Number of ROIs
    tf = 0;
    atlas_choices = {
        'AAL3'
        'Schaeffer 400 (17 networks)'
        'Schaeffer 200 (17 networks)'
        'Schaeffer_neurosynth'
        };
    atlas_prompt = 'Which atlas was used in CONN preprocessing?';
    while tf==0
        if ~isfield(dFC_project,'atlas')
            [atlas,tf] = listdlg('PromptString',atlas_prompt,'SelectionMode','single','ListString',atlas_choices);
            if tf
                dFC_project.atlas = atlas;
            end
        end
    end
    switch dFC_project.atlas
        case 1
            tf = 0;
            network_choices = {
                'Language Network'
                'Salience Network'
                'Dorsal Attention Network'
                'Default Mode Network'
                'Random Network'
                'Whole Brain'
                'Whole Cortex'
                'Select individual ROIs'
                };
            network_prompt = 'Which ROIs would you like to include? (AAL3 Atlas)';
            while tf==0
                if ~isfield(dFC_project,'include')
                    [network,tf] = listdlg('PromptString',network_prompt,'SelectionMode','single','ListString',network_choices);
                    if tf==0 || network==6
                        dFC_project.include = [1:nroi];
                    elseif network == numel(network_choices)
                        dFC_project = get_roi_info(dFC_project);
                        [dFC_project.include,tf] = listdlg('ListString',dFC_project.all_roi_names);
                    else
                        dFC_project = get_network(dFC_project,network);
                    end
                    if dFC_project.which_hemi==2
                        left = find( ...
                            [1	0	1	0	1	0	1	0	1	0	1	0 ...
                            1	0	1	0	1	0	1	0	1	0	1	0 ...
                            1	0	1	0	1	0	1	0	1	0	1	0 ...
                            1	0	1	0	1	0	1	0	1	0	1	0 ...
                            1	0	1	0	1	0	1	0	1	0	1	0 ...
                            1	0	1	0	1	0	1	0	1	0	1	0 ...
                            1	0	1	0	1	0	1	0	1	0	1	0 ...
                            1	0	1	0	1	0	1	0	1	0	1	0 ...
                            1	0	1	0	1	0	1	0	1	0	1	0 ...
                            0	0	0	0	0	0	0	0	1	0	1	0 ...
                            1	0	1	0	1	0	1	0	1	0	1	0 ...
                            1	0	1	0	1	0	1	0	1	0	1	0 ...
                            1	0	1	0	1	0	1	0	1	0	1	0 ...
                            1	0	1	0	1	0	1	0	0	0]);
                        dFC_project.include = dFC_project.include(ismember(dFC_project.include,left));
                        dFC_project.which_hemi = 1;
                    elseif dFC_project.which_hemi==3
                        right = find( ...
                            [0	1	0	1	0	1	0	1	0	1	0	1 ...
                            0	1	0	1	0	1	0	1	0	1	0	1 ...
                            0	1	0	1	0	1	0	1	0	1	0	1 ...
                            0	1	0	1	0	1	0	1	0	1	0	1 ...
                            0	1	0	1	0	1	0	1	0	1	0	1 ...
                            0	1	0	1	0	1	0	1	0	1	0	1 ...
                            0	1	0	1	0	1	0	1	0	1	0	1 ...
                            0	1	0	1	0	1	0	1	0	1	0	1 ...
                            0	1	0	1	0	1	0	1	0	1	0	1 ...
                            0	0	0	0	0	0	0	0	0	1	0	1 ...
                            0	1	0	1	0	1	0	1	0	1	0	1 ...
                            0	1	0	1	0	1	0	1	0	1	0	1 ...
                            0	1	0	1	0	1	0	1	0	1	0	1 ...
                            0	1	0	1	0	1	0	1	0	0]);
                        dFC_project.include = dFC_project.include(ismember(dFC_project.include,right));
                        dFC_project.which_hemi = 1;
                    end
                end
            end
        case 2
            network_choices = get_schaefer_networks(dFC_project);
            network_choices{numel(network_choices)+1} = 'All ROIs';
            network_prompt = 'Which network(s) would you like to analyze?';
            tf = 0;
            while tf==0
                if ~isfield(dFC_project,'include')
                    [network,tf] = listdlg('PromptString',network_prompt,'SelectionMode','multiple','ListString',network_choices);
                    if tf==0
                        msgbox('You must make a selection')
                    elseif numel(network)==1 && network == numel(network_choices)
                        dFC_project.include = [1:nroi];
                        dFC_project = get_roi_info(dFC_project);
                    elseif numel(network)>1 && ismember(numel(network_choices),network)
                        msgbox('"All ROIs" cannot be combined with other network choices')
                        tf = 0;
                    else
                        dFC_project = get_network(dFC_project,network);
                        dFC_project = get_roi_info(dFC_project);
                    end
                end
            end
        case 3
            tf = 0;
            network_choices = get_schaefer_networks(dFC_project);
            network_choices{numel(network_choices)+1} = 'All ROIs';
            network_prompt = 'Which network(s) would you like to analyze?';
            while tf==0
                if ~isfield(dFC_project,'include')
                    [network,tf] = listdlg('PromptString',network_prompt,'SelectionMode','multiple','ListString',network_choices);
                    if tf==0
                        msgbox('You must make a selection')
                    elseif numel(network)==1 && network == numel(network_choices)
                        dFC_project.include = [1:nroi];
                        dFC_project = get_roi_info(dFC_project);
                    elseif numel(network)>1 && ismember(numel(network_choices),network)
                        msgbox('"All ROIs" cannot be combined with other network choices')
                        tf = 0;
                    else
                        dFC_project = get_network(dFC_project,network);
                    end
                end
            end
        case 4
            network_choices = {
                'Language Network'
                'Salience Network'
                'Dorsal Attention Network'
                'Default Mode Network'
                'Select individual ROIs'
                };
            network_prompt = 'Which network/ROIs do you want to include?';
            tf = 0;
            while tf==0
                if ~isfield(dFC_project,'include')
                    [network,tf] = listdlg('SelectionMode','single','ListString',network_choices);
                    if tf==0
                        return
                    elseif network == numel(network_choices)
                        dFC_project = get_roi_info(dFC_project);
                        [dFC_project.include,tf] = listdlg('PromptString',network_prompt,'ListString',dFC_project.all_roi_names);
                    else
                        dFC_project = get_network(dFC_project,network);
                    end
                end
                dFC_project = get_roi_info(dFC_project);
            end
        otherwise
            disp oops!
    end
        
    dFC_project = get_roi_info(dFC_project);
    disp('Getting number of volumes for each subject...')
    for i=1:nsub
        for k=1:ncon_actual
            volik = load([sprintf('ROI_Subject%03d_',i),sprintf('Condition%03d.mat',dFC_project.conds(k))]).data{full_roi-ncon+dFC_project.conds(k)} ~= 0;
            nvols(i,k) = length(volik);
            nvols_actual(i,k) = length(volik(volik~=0));
%             disp(['Number of volumes = ',num2str(nvols_actual(i,k))])
        end
    end
    nvols_lin = reshape(nvols,nsub*ncon_actual,1);
    nvols_actual_lin = reshape(nvols_actual,nsub*ncon_actual,1);
    nvol = max(nvols_lin);
    nvol_cut = min(nvols_actual_lin);
    disp(['Minimum number of volumes is ',num2str(nvol_cut),'. ',num2str(nvol_cut),' volumes for each subject will be included in the analysis'])
    ROIs = load('ROI_Subject001_Condition001.mat').names(4:(nroi+3));
    ROIs = ROIs(dFC_project.include);
end

disp(['Extracting time series data for ',num2str(numel(dFC_project.include)),' ROIs'])

if isfield(dFC_project,'dur')
    nvol = dFC_project.dur*60/min(dFC_project.TR);
end
BOLD = zeros(nvol,numel(dFC_project.include),nsub,ncon_actual);


if numel(dFC_project.TR) == 1
    for k=1:ncon_actual
        files = dir(sprintf('*_Condition%03d.mat',dFC_project.conds(k)));
        
        for i=1:numel(files)
        
            clear D
            data = load(files(i).name).data;
            con_check = load(files(i).name).data{full_roi-ncon+dFC_project.conds(k)} ~= 0;
            
            D = zeros(length(data{1}),length(dFC_project.include));
            
            for j=1:length(dFC_project.include)
                D(:,j) = data{dFC_project.include(j)+3};
            end
            
    %         BOLD(:,:,i,k) = D;
    %         disp(['subject = ',num2str(i)])
            BOLD(1:sum(con_check),:,i,k) = D(con_check,:);
            disp(['Time series data for Subject ',num2str(i),' Condition ',num2str(dFC_project.conds(k)),' extracted ...'])
    
        end
    end

    BOLD = BOLD(1:nvol_cut,:,:,:);
else
    dFC_project.all_nvols = zeros(nsub,1);
    for k=1:ncon_actual
        files = dir(sprintf('*_Condition%03d.mat',dFC_project.conds(k)));
        
        for i=1:numel(files)
        
            clear D
            data = load(files(i).name).data;
            con_check = load(files(i).name).data{full_roi-ncon+dFC_project.conds(k)} ~= 0;
            
            D = zeros(length(data{1}),length(dFC_project.include));
            
            for j=1:length(dFC_project.include)
                D(:,j) = data{dFC_project.include(j)+3};
            end

            if ~ismember(max(dFC_project.TR)*sum(con_check)/nvol,dFC_project.TR)
                disp(['TR for current subject data appears to be ',num2str(max(dFC_project.TR)*length(data{1})/nvol),' which does not match any provided TR values'])
            end
            
    %         BOLD(:,:,i,k) = D;
    %         disp(['subject = ',num2str(i)])
            if k==1
                dFC_project.all_nvols(i) = sum(con_check);
            end
            D = D(con_check,:);
            D = imresize(D,[nvol,size(D,2)]);
            BOLD(:,:,i,k) = D;
            disp(['Time series data for Subject ',num2str(i),' Condition ',num2str(dFC_project.conds(k)),' extracted ...'])
    
        end
    end
    BOLD = BOLD(1:nvol,:,:,:);
end

% Sort BOLD by hemisphere
[dFC_project.hemispheres,sort_index] = sort(dFC_project.hemispheres);
dFC_project.roi_names = dFC_project.roi_names(sort_index);
BOLD = BOLD(:,sort_index,:,:);
ROIs = ROIs(sort_index);
dFC_project.include_hemi_sorted = dFC_project.include(sort_index);

if dFC_project.which_hemi == 2
    retain = cellfun(@(x)strcmp(x,'L'),dFC_project.hemispheres);
elseif dFC_project.which_hemi == 3
    retain = cellfun(@(x)strcmp(x,'R'),dFC_project.hemispheres);
end

if dFC_project.which_hemi ~= 1
    dFC_project.roi_names = dFC_project.roi_names(retain);
    BOLD = BOLD(:,retain,:,:);
    ROIs = ROIs(retain);
    dFC_project.hemispheres = dFC_project.hemispheres(retain);
end

if numel(dFC_project.TR) == 1
    nvol = nvol_cut;
end



cd(oldPath)

end
