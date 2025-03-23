function FC_project = get_roi_info(FC_project)

old_path = pwd;

switch FC_project.atlas
    case 1 % AAL3 parcellation
        roi_info = struct2cell(load('ROI_MNI_V7_List.mat').ROI);
        all_roi_names = permute(roi_info(3,1,:),[3 1 2]);
        
        % Get hemisphere info
        roi_split = cellfun(@(x)split(x,'_'),all_roi_names,'UniformOutput',false); % Split ROI names at '_'
        all_hemispheres = cellfun(@(v) v(end),roi_split); % Get hemisphere from last element of ROI name
    case 2 % Schaefer 400 parcellation
        cd(FC_project.CONN_results_dir)
        all_roi_names = schaefer_roi_names();
        all_hemispheres = cell(numel(all_roi_names),1);
        for i=1:numel(all_roi_names)
            split_name = split(all_roi_names{i},'_');
            hemi_name = split_name{2};
            all_hemispheres{i} = hemi_name(1);
        end
    case 3 % Schaefer 200 parcellation
        cd(FC_project.CONN_results_dir)
        all_roi_names = schaefer_roi_names_200();
        all_hemispheres = cell(numel(all_roi_names),1);
        for i=1:numel(all_roi_names)
            split_name = split(all_roi_names{i},'_');
            hemi_name = split_name{2};
            all_hemispheres{i} = hemi_name(1);
        end
    case 4
        path = pwd;
        cd('C:\Users\isaac\OneDrive\Documents\Grad School\Kiran Lab\DIssertation\Dynamic Functional Connectivity Custom MATLAB Scripts\schaefer_neurosynth_ROIs')
        files = dir('CBR_neurosynth_*.nii');
        all_roi_names = transpose({files.name});
        roi_split = cellfun(@(x)split(x,'_'),all_roi_names,'UniformOutput',false); % Split ROI names at '_'
        all_hemispheres = cellfun(@(v) v(end),roi_split); % Get hemisphere from last element of ROI name
        roi_split = cellfun(@(x)split(x,'.'),all_hemispheres,'UniformOutput',false); % Split ROI names at '_'
        all_hemispheres = cellfun(@(v) v(1),roi_split); % Get hemisphere from last element of ROI name
        cd(path)
    otherwise
        disp oops!
end

for i=1:numel(all_hemispheres)
    if ~strcmp(all_hemispheres{i},'L') && ~strcmp(all_hemispheres{i},'R')
        all_hemispheres{i} = 'M';
    end
end

if isfield(FC_project,'include')
    FC_project.hemispheres = all_hemispheres(FC_project.include);
    FC_project.roi_names = all_roi_names(FC_project.include);
end

FC_project.all_hemispheres = all_hemispheres;
FC_project.all_roi_names = all_roi_names;

cd(old_path)

end