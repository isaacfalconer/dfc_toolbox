dFC_project = struct();
[dFC_project] = compute_dFC(dFC_project,1);
if ~isstruct(dFC_project)
    errordlg('An error occured while computing dFC')
    return
end
dFC_project = compute_variability(dFC_project);

if ~isfield(dFC_project,'var_output_dir')
    dFC_project.var_output_dir =uigetdir(pwd,'Select folder for writing results');
end
old_path = pwd;
cd(dFC_project.var_output_dir)
writematrix(dFC_project.mean_var,'mean_var.csv')
writematrix(dFC_project.TF_comm_struct,'TF_comm_struct.csv')
writematrix(dFC_project.mean_similarity,'mean_similarity.csv')

if isfield(dFC_project,'proportion_spared')
    proportion_spared = dFC_project.proportion_spared(dFC_project.include_hemi_sorted,:);
    voxels_spared = dFC_project.voxels_spared(dFC_project.include_hemi_sorted,:);
    total_voxels_spared = transpose(sum(voxels_spared,1));
    writematrix(proportion_spared,'proportion_spared.csv')
    writematrix(total_voxels_spared,'total_voxels_spared.csv')
    writematrix(dFC_project.lesion_volume,'lesion_volume.csv')
end


% writematrix(dFC_project.mean_duration,'mean_duration.csv')
% writematrix(dFC_project.time_in_state,'time_in_state.csv')
% writematrix(dFC_project.mean_in_state_dFC,'mean_in_state_dFC.csv')
% writematrix(dFC_project.mean_out_of_state_dFC,'mean_out_of_state_dFC.csv')
for i=1:numel(dFC_project.vars)
    writematrix(dFC_project.vars{i},['vars_sub_',num2str(i),'.csv'])
end

dFC_project = get_roi_lobe(dFC_project);
writecell(dFC_project.ROI_location,'ROI_location.csv')
writecell(dFC_project.hemispheres,'hemispheres.csv')

if ~isfield(dFC_project,'project_name')
    dFC_project.project_name = char(inputdlg('Enter a name for your project ', 'Save Project...', [1 45]));
end
save_path = which('findme');
cd(save_path(1:length(save_path)-9));

if exist([dFC_project.project_name,'.mat'],'file')
    delete([dFC_project.project_name,'.mat'])
end

x = app.center_x;
y = app.center_y;
fig = uifigure('Position',[x-85 y-35 170 70]);
uitextarea(fig,'Value','Saving Project. Please wait...','Position',[10 10 150 50]);
drawnow
figure(app.DynamicFunctionalConnectivityToolboxUIFigure)
figure(fig)
pause(0.01)

if whos('dFC_project').bytes<2000000000
    save([dFC_project.project_name,'.mat'], 'dFC_project')
else
    save([dFC_project.project_name,'.mat'], 'dFC_project','-v7.3')
end

close(fig)

cd(old_path)