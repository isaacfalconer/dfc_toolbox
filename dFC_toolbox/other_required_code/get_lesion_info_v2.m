function dFC_project = get_lesion_info_v2(dFC_project)

old_path = pwd;

cd(fileparts(which('AAL3v1.nii')))
AAL3 = load_nii('AAL3v1.nii');
AAL3_origin = AAL3.hdr.hist.originator(1:3);
AAL3_img = double(AAL3.img);
AAL3_img(isnan(AAL3_img)) = 0;
cd('C:\Users\isaac\OneDrive\Documents\Grad School\Kiran Lab\DIssertation\Aim 2\Data\lesion_maps')
check = 0;
while check==0
    [lesion_files,lesion_path] = uigetfile('*.nii*',  'All Files (*.*)','MultiSelect','on');
    if ~ischar(lesion_files) && ~iscell(lesion_files) && lesion_files==0
        return
    elseif numel(lesion_files) ~= dFC_project.nsub
        waitfor(warndlg('The number of lesion map files must match the number of subjects'))
    else
        check = 1;
    end
end
threshold = inputdlg('What threshold would you like to use to exclude ROIs? (proportion damaged)','',[1 40],{'0.5'});
threshold = eval(threshold{1});
cd(lesion_path)

% lesions = zeros([size(AAL3_img) numel(lesion_files)]);
proportion_spared = zeros(166,numel(lesion_files));
voxels_spared = zeros(166,numel(lesion_files));
excluded_ROIs = zeros(166,numel(lesion_files));
lesion_volume = zeros(numel(lesion_files),1);
for i=1:numel(lesion_files)
    disp(['Importing lesion data for subject ',num2str(i),'...'])
    lesion = load_nii(lesion_files{i});
    if any(lesion.hdr.dime.pixdim(2:4)~=[2 2 2])
        disp(['Resampling subject ',num2str(i),' lesion map to 2x2x2mm voxel size and writing to new file...'])
        new_name = ['2mm_',lesion_files{i}];
        reslice_nii(lesion_files{i},new_name,[2 2 2]);
        lesion = load_nii(new_name);
    end
    lesion_origin = lesion.hdr.hist.originator(1:3);
    origin_diffs = AAL3_origin - lesion_origin;
    if origin_diffs(1)>=0 
        opt.pad_from_L = origin_diffs(1);
        opt.pad_from_P = origin_diffs(2);
        opt.pad_from_I = origin_diffs(3);
        lesion = pad_nii(lesion, opt);
        lesions(:,:,:,i) = lesion.img;
        lesions(isnan(lesions)) = 0;
        for j=1:166
            roi_lesion_overlap = (AAL3_img==convert_2_AAL3_indices(j)) .* lesions(:,:,:,i);
            pct_damage = sum(roi_lesion_overlap,'all')/sum(AAL3_img==convert_2_AAL3_indices(j),'all');
            if pct_damage>=threshold
                excluded_ROIs(j,i) = 1;
            end
            proportion_spared(j,i) = 1-pct_damage;
            voxels_spared(j,i) = sum(AAL3_img==convert_2_AAL3_indices(j),'all') - sum(roi_lesion_overlap,'all');
        end
        lesion_volume(i) = sum(lesions(:,:,:,i),'all');
        % old_path2 = pwd;
        % cd('C:\Users\isaac\OneDrive\Documents\Grad School\Kiran Lab\DIssertation\Aim 1\Code\lesion_code_testing')
        % roi_labeled_lesion = AAL3;
        % roi_labeled_lesion(lesions(:,:,:,i)~=1) = 0;
        % niftiwrite(flip(roi_labeled_lesion,1),['roi_labeled_',lesion_files{i}])
        % cd(old_path2)
    elseif origin_diffs(1)<0
        origin_diffs = lesion_origin - AAL3_origin;
        opt.pad_from_L = origin_diffs(1);
        opt.pad_from_P = origin_diffs(2);
        opt.pad_from_I = origin_diffs(3);
        AAL3_temp = pad_nii(AAL3, opt);
        lesions(:,:,:,i) = lesion.img;
        lesions(isnan(lesions)) = 0;
        for j=1:166
            x = 1:size(AAL3_temp.img,1);
            y = 1:size(AAL3_temp.img,2);
            z = 1:size(AAL3_temp.img,3);
            roi_lesion_overlap = (AAL3_temp.img==convert_2_AAL3_indices(j)) .* lesions(x,y,z,i);
            pct_damage = sum(roi_lesion_overlap,'all')/sum(AAL3_temp.img==convert_2_AAL3_indices(j),'all');
            if pct_damage>=threshold
                excluded_ROIs(j,i) = 1;
            end
            proportion_spared(j,i) = 1-pct_damage;
            voxels_spared(j,i) = sum(AAL3_temp.img==convert_2_AAL3_indices(j),'all') - sum(roi_lesion_overlap,'all');
        end
        lesion_volume(i) = sum(lesions(:,:,:,i),'all');
        % old_path2 = pwd;
        % cd('C:\Users\isaac\OneDrive\Documents\Grad School\Kiran Lab\DIssertation\Aim 1\Code\lesion_code_testing')
        % roi_labeled_lesion = AAL3;
        % roi_labeled_lesion(lesions(:,:,:,i)~=1) = 0;
        % niftiwrite(flip(roi_labeled_lesion,1),['roi_labeled_',lesion_files{i}])
        % cd(old_path2)
    else
        warndlg("Lesion import failed")
    end
end

cd(old_path);

dFC_project.lesions = lesions;
dFC_project.excluded_ROIs = excluded_ROIs;
dFC_project.proportion_spared = proportion_spared;
dFC_project.voxels_spared = voxels_spared;
dFC_project.lesion_volume = lesion_volume;

end