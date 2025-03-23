function dFC_project = compute_dFC(dFC_project,new_project)


if ~isfield(dFC_project,'CONN_preprocessing_dir')
    dFC_project.CONN_preprocessing_dir = uigetdir(pwd,'Select folder containing output of CONN preprocessing');
end
cd(dFC_project.CONN_preprocessing_dir)
cd ..

if ~isfield(dFC_project,'CONN_results_dir')
    dFC_project.CONN_results_dir = uigetdir(pwd,'Select folder containing output of CONN first level analyses');
end

% old_path = pwd;
% cd(dFC_project.CONN_results_dir)
% for i=1:

dFC_project = get_ncon_nsub(dFC_project);
dFC_project = get_nroi(dFC_project);

if new_project
    one_TR = questdlg('Do all scans have the same TR?');
    one_dur = questdlg('Do all scans have the same duration?');
    if strcmp('Yes',one_TR)
        prompt = {'Enter desired window length in seconds:'
            'Enter repetition time of functional scan (TR):'
            };
        dlgtitle = 'Input';
        dims = [1 75];
        defInput = {'30','2.4'};
        responses = inputdlg(prompt,dlgtitle,dims,defInput);
        if isempty(responses)
            dFC_project = 0;
            return
        end
        dFC_project.window = str2double(responses{1});
        dFC_project.TR = str2double(responses{2});
    elseif strcmp('No',one_TR)
        TR_str = inputdlg('Enter all unique TR values separated by commas (e.g., 2,2.4):','Repitition Times');
        dFC_project.TR = str2double(split(TR_str,','));
        if strcmp('No',one_dur)
            disp Durations must be consistent if TR values vary.
            return
            % dFC_project.dur = zeros(1,length(dFC_project.TR));
            % for i=1:length(dFC_project.TR)
            %     dur_str = inputdlg(['Enter duration (in minutes) for scans with TR of ',num2str(dFC_project.TR(i))],'Scan Duration');
            %     dFC_project.dur(i) = str2double(split(dur_str,','));
            % end
        else
            dur_str = inputdlg('Enter scan duration (in minutes):','Scan Duration');
            dFC_project.dur = str2double(split(dur_str,','));
        end
        window_str = inputdlg('Enter desired window length in seconds','Scan Durations');
        dFC_project.window = str2double(split(window_str,','));
    else
        disp oops!
    end
end

dFC_project = get_conds(dFC_project);

if ~isfield(dFC_project,'cond_names')
    dFC_project.cond_names = cell(numel(dFC_project.conds),1);
    for i=1:numel(dFC_project.conds)
        dFC_project.cond_names{i} = char(inputdlg(['Enter name of Condition ',num2str(dFC_project.conds(i))], ""));
    end
end


dFC_project.which_hemi = menu('Which hemisphere(s) would you like to include?','Left and Right','Left Only','Right Only');

window = floor(dFC_project.window/min(dFC_project.TR)); % convert window length from seconds to number of volumes

[BOLD,ROIs,nvol,dFC_project] = read_BOLD(dFC_project);
dFC_project.BOLD = BOLD;
nroi = numel(ROIs);
nwin = nvol-(window-1);
ncon = numel(dFC_project.cond_names);
nsub = dFC_project.nsub;
% ncol = ((nroi^2)-nroi)/2-(nroi/2)^2+nroi/2; % number of intrahemisphereic and homotopic connections

if dFC_project.atlas==1
    lesion_prompt = ['Do you have lesion map files to import? (Damaged ROIs will' ...
        ' be excluded from certain analyses using a threshold of your choice)'];
    if strcmp(questdlg(lesion_prompt),'Yes')
        dFC_project = get_lesion_info_v2(dFC_project);
    else
        dFC_project.excluded_ROIs = zeros(numel(dFC_project.all_roi_names),nsub);
    end
end

fig = uifigure;
progress_bar = uiprogressdlg(fig,'Title','Computing dFC. This may take several minutes...','Message','');
pause(0.01)

dFC = zeros(nroi,nroi,nwin,ncon,nsub);
% dFC_flat = zeros(nwin,ncol,nsub);
% columns = cell(ncol,1);
% col = 1;

% Compute dFC
for i=1:nsub
    for j=1:ncon
        % disp(['Computing dFC for Subject ',num2str(i),' Condition ',num2str(j),'...'])
        progress_bar.Message = ['Computing dFC for Subject ',num2str(i),' Condition ',num2str(j),'...'];
        progress_bar.Value = (i-1)/nsub + (j/ncon)*(1/nsub);
        pause(0.01)
        for k=1:nwin
            for n=1:nroi
                for m=1:nroi
                    if m>n % && (strcmp(dFC_project.hemispheres{m},dFC_project.hemispheres{n}) || strcmp(dFC_project.roi_names{m},dFC_project.roi_names{n}))
                        corr = corrcoef(BOLD(k:(k+window-1),m,i,j),BOLD(k:(k+window-1),n,i,j));
                        dFC(m,n,k,j,i) = corr(2,1);
                        dFC(n,m,k,j,i) = corr(2,1);
                    end
                end
            end
            col = 1;
        end
    end
end

dFC(isnan(dFC)) = 0;

% % Binarize connectivity
% dFC(dFC >= 0.4) = 1;
% dFC(dFC < 0.4) = 0;

dFC = atanh(dFC);

dFC(isnan(dFC)) = 0;

dFC_project.dFC = dFC;
% dFC_project.columns = columns;
% dFC_project.dFC_flat = dFC_flat;

% dFC_project = intersub_var(BOLD,dFC_project);
% 
% dFC_project.dFC(dFC_project.intersub_var >= median(dFC_project.intersub_var)) = 0;

% dFC_project = compute_TF_FCD(dFC_project,1);

close(progress_bar)
close(fig)

dFC_project = compute_TF_comm_struct_v4(dFC_project,1);

end