function [dFC_project] = dFC_kmeans_app(dFC_project, condition, maxk)

rng_state = rng();


% % Ensure dFC matrix has 5 dimensions
% if numel(d)~=5 && dFC_project.nsub~=1
%     disp 'dFC matrix has wrong number of dimensions'
%     return
% end

if strcmp(questdlg('Do you want limit k-means clustering to right hemisphere data?'),'Yes')
    retain = cellfun(@(x)strcmp(x,'R'),dFC_project.hemispheres);
    dFC = dFC_project.dFC(retain,retain,:,condition,:); % Only use specified condition
else
    dFC = dFC_project.dFC(:,:,:,condition,:); % Only use specified condition
end

% dFC = dFC_project.dFC_comm_struct;
% mean_dFC = mean(dFC(dFC>0));
% dFC(dFC==0) = mean_dFC;
% for i=1:size(dFC,1)
%     dFC(i,i,:,:) = 0;
% end
% dFC_project.dFC_states_analyses = dFC;
choice = 'No';
if strcmp(choice,'Yes')
    d = size(dFC_project.dFC_for_kmeans); 
    nwin = d(3);                % Get number of volumes in dFC matrix
    nroi = d(1);    % Get number of rois
    nsub = dFC_project.nsub;    % Get number of subjects
    dFC_flat = dFC_project.dFC_flat;
elseif strcmp(choice,'No')
    d = size(dFC); 
    nwin = d(3);                % Get number of volumes in dFC matrix
    nroi = d(1);    % Get number of rois
    nsub = dFC_project.nsub;    % Get number of subjects
    dFC_perm = permute(dFC,[3 5 1 2 4]);             % Permute matrix to eliminate condition dimension
    % dFC_perm = permute(dFC,[3 4 1 2]);
    dFC_flat = reshape(dFC_perm,[size(dFC,3)*nsub nroi^2]); % Reshape matrix to an NxM matrix with earch row corresponding to one volume and one subjectand columns corresponding to ROI pairs
    else
    disp oops!
end
% dFC_flat(dFC_flat<clust_threshold) = 0;


ks = (transpose([1,2,3,4,6,8,10,12,20,28,36,44,52,68,84,100]));
ks = [ks(ks<maxk);maxk];

IDX = zeros(numel(dFC_flat(:,1)),numel(ks));
C = cell(numel(ks),1);
swss = zeros(numel(ks),1);

rng(123)

if strcmp(questdlg('Do you want to perform a new PCA?'),'Yes')
    [~,pca_score,pca_latent] = pca(dFC_flat);
    if pca_latent(2)>=1
        keep = pca_latent>=1;
        if length(keep)>50
            keep = 1:50;
        end
    else
        keep=[1,2];
    end
    dFC_pca = pca_score(:,keep);
else
    dFC_pca = dFC_project.dFC_pca;
    pca_score = dFC_pca;
end

dFC_project.dFC_pca = dFC_pca;
if nroi>50
    clust_data = dFC_pca;
else
    clust_data = dFC_flat;
end


for i=1:numel(ks)
    rng(123)
    [IDX(:,i), C{i}, sumD] = kmeans(clust_data,ks(i));
    swss(i,1) = sum(sumD);
    disp(['Clustering for k = ',num2str(ks(i)),' done...'])
end

[xData, yData] = prepareCurveData( ks, swss );
dFC_project.elbow_plot_k = xData;
dFC_project.elbow_plot_swss = yData;

% Plot fit with data.
fig = figure( 'Name', 'Elbow Plot' );
% h = plot( fitresult, xData, yData );
plot(xData,yData)
% legend( h, 'Total WSS vs. k', 'Elbow Plot', 'Location', 'NorthEast', 'Interpreter', 'none' );
% Label axes
xlabel( 'Number of Clusters', 'Interpreter', 'none' );
ylabel( 'Total Within-Cluster Sum of Squares', 'Interpreter', 'none' );
grid on


% disp(gof)

prompt = {'Minimum k:','Maximum k:'};
dlgtitle = 'Range of k values for silhouette analysis';
dims = [1 85];
k_range = inputdlg(prompt,dlgtitle,dims);
numk = (k_range{2}-k_range{1}+1);
k_min = str2double(k_range{1});
k_max = str2double(k_range{2});
k_sil = [k_min:k_max];

close(fig)

fig = figure('Position',[700,200,200*ceil(numel(k_sil)/floor(sqrt(numel(k_sil)))),200*floor(sqrt(numel(k_sil)))]);
tiledlayout(floor(sqrt(numel(k_sil))),ceil(numel(k_sil)/floor(sqrt(numel(k_sil)))))
for i=1:numel(k_sil)

    rng(123)
    clust = kmeans(clust_data,k_sil(i));
    % colors = {'black','red','green','blue','magenta','yellow','cyan',[0.8500 0.3250 0.0980],[0.4940 0.1840 0.5560],[0.6350 0.0780 0.1840]};
    colors = hsv;
    nexttile
    hold on
    for j=1:k_sil(i)
        % DT = DelaunayTri(pca_score(clust==j,1:3));
        % hullFacets = convexHull(DT);
        % trisurf(hullFacets,DT.X(:,1),DT.X(:,2),DT.X(:,3),'FaceColor',colors(j*floor(256/k_sil(i)),:),'FaceAlpha',0.2,'EdgeAlpha',0.2)
        % view(0,90)
        % scatter3(pca_score(clust==j,1),pca_score(clust==j,2),pca_score(clust==j,3),1,colors(j*floor(256/k_sil(i)),:),'filled')
        scatter(dFC_pca(clust==j,1),dFC_pca(clust==j,2),1,colors(j*floor(256/k_sil(i)),:),'filled')
    end
    hold off
    title(['PCA plot for k = ',num2str(k_sil(i))])
    xlabel('PC1')
    ylabel('PC2')

end




k_opt = eval(cell2mat(inputdlg('Enter desired number of clusters according to PCA plots')));

close(fig)

rng(123)
[dFC_project.clusters , ~, ~] = kmeans(clust_data,k_opt);
for i=1:nsub

end

state_dFC = cell(k_opt,1);
dFC_project.state_dFC_mean = cell(k_opt,1);

% Get directory containing BOLD time series .mat files
if ~isfield(dFC_project,'CONN_results_dir')
    dFC_project.CONN_results_dir = uigetdir(pwd,'Select folder containing output of CONN first level analyses');
end

for i=1:k_opt
    state_dFC{i} = [];
    if strcmp(choice,'Yes')
        dFC_project.state_dFC_mean{i} = mean(dFC_project.dFC_for_kmeans(:,:,dFC_project.clusters==i),[3]);
    elseif strcmp(choice,'No')
        try
            full_dFC = dFC_project.dFC;
            full_nroi = size(full_dFC,1);
            dFC_cat = reshape(full_dFC,[full_nroi full_nroi size(dFC,3)*nsub]);
        catch
            waitfor(errordlg('dFC_cat is the problem'))
        end
        dFC_project.state_dFC_mean{i} = mean(dFC_cat(:,:,dFC_project.clusters==i),[3]);
    end
end

% dFC_project.states_mean = mean(cat(3,dFC_project.state_dFC_mean{:}),3);

dFC_project.PC1 = dFC_pca(:,1);
dFC_project.PC2 = dFC_pca(:,2);
dFC_project.state_occ = zeros(nsub,k_opt);

if strcmp(choice,'Yes')
    for i=1:k_opt
        for j=1:nsub
            sub_clust = dFC_project.clusters(dFC_project.subs_for_kmeans==j);
            dFC_project.state_occ(j,i) = sum(sub_clust==i);
        end
    end
elseif strcmp(choice,'No')
    disp('179')
    sub_clust = reshape(dFC_project.clusters,[size(dFC,3) nsub]);
    for i=1:k_opt
        for j=1:nsub
            dFC_project.state_occ(j,i) = sum(sub_clust(:,j)==i);
        end
    end
else
    disp oops!
end


dFC_project.TF_gr = compute_TF(sub_clust, min(dFC_project.TR));

path = pwd;
dFC_project.var_output_dir = uigetdir(pwd,'Select folder for writing results');
cd(dFC_project.var_output_dir)
writematrix(dFC_project.state_occ,'state_occ.csv')
writematrix(dFC_project.TF_gr,'TF_gr.csv')
cd(path)

rng(rng_state)

end