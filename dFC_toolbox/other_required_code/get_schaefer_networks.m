function network_names = get_schaefer_networks(FC_project)

old_path = pwd;

cd(FC_project.CONN_results_dir)
names = schaefer_roi_names();
network_names = cell(numel(names),1);
for i=1:numel(names)
    split_name = split(names{i},'_');
    network_names{i} = split_name{3};
end

network_names = unique(network_names);
cd(old_path)

end