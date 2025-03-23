legend(app.UIAxes,'off')

nroi = size(app.dFC_project.dFC,1);

default_colors = transpose([(257-[1:256])/256; [1:128]/128,(129-[1:128])/128; [1:256]/256]);
default = [0 0.4470 0.7410;
0.900 0.150 0.0980;
0.9290 0.6940 0.1250;
0.4940 0.1840 0.5560;
0.2660 0.6740 0.4880;
0.3010 0.7450 0.9330;
0.6350 0.0780 0.1840;
0 0.4470 0.7410;
0.8500 0.3250 0.0980;
0.9290 0.6940 0.1250;
0.4940 0.1840 0.5560;
0.4660 0.6740 0.1880;
0.3010 0.7450 0.9330;
0.6350 0.0780 0.1840;
0 0.4470 0.7410;
0.8500 0.3250 0.0980;
0.9290 0.6940 0.1250;
0.4940 0.1840 0.5560;
0.4660 0.6740 0.1880;
0.3010 0.7450 0.9330;
0.6350 0.0780 0.1840;
0 0.4470 0.7410;
0.8500 0.3250 0.0980;
0.9290 0.6940 0.1250;
0.4940 0.1840 0.5560;
0.4660 0.6740 0.1880;
0.3010 0.7450 0.9330;
0.6350 0.0780 0.1840
0 0.4470 0.7410;
0.8500 0.3250 0.0980;
0.9290 0.6940 0.1250;
0.4940 0.1840 0.5560;
0.4660 0.6740 0.1880;
0.3010 0.7450 0.9330;
0.6350 0.0780 0.1840;
0 0.4470 0.7410;
0.8500 0.3250 0.0980;
0.9290 0.6940 0.1250;
0.4940 0.1840 0.5560;
0.4660 0.6740 0.1880;
0.3010 0.7450 0.9330;
0.6350 0.0780 0.1840];

rycb = flip(transpose([ones(1,85) (86-[1:85])/85 zeros(1,86);
                 [1:85]/85 ones(1,85) (87-[1:86])/86;
                 zeros(1,86) [1:85]/85 ones(1,85)]),1);

redblue = transpose([(257-[1:256])/256;
                     0.2*ones(1,256);
                     [1:256]/256]);

earth = transpose([(129-[1:128])/151+19/128,(19/128)*ones(1,128); ...
                           [1:128]/210+46/256,(129-[1:128])/210+46/256; ...
                           (19/128)*ones(1,128),[1:128]/151+19/128]);

if app.dFC_project.atlas==1 && (nroi==166 || (nroi == 78 && app.dFC_project.which_hemi>1))
    % rois = 1:166;
    rois = app.dFC_project.include_hemi_sorted(ismember(app.dFC_project.include_hemi_sorted,rois));
    % if numel(rois)==166 && app.dFC_project.which_hemi==2
    %     rois = rois(1:78);
    %     hemi = 1:78;
    % elseif numel(rois)==166 && app.dFC_project.which_hemi==3
    %     rois = rois(89:166);
    %     hemi = 89:166;
    % % elseif numel(rois)==166
    % %     hemi = 1:166;
    % elseif app.dFC_project.which_hemi==2
    %     rois = rois(1:(numel(rois)/2));
    %     hemi = 1:78;
    % elseif app.dFC_project.which_hemi==3
    %     rois = rois((numel(rois)/2) + (1:(numel(rois)/2)));
    %     hemi = 89:166;
    % else
    %     hemi = 1:166;
    % end

    mygraph = app.dFC_project.state_dFC_mean{state};
    % app.dFC_project.include_hemi_sorted(hemi)
    % size(mygraph)
    % size(ismember(app.dFC_project.include_hemi_sorted(hemi),rois))
    mygraph = mygraph(ismember(app.dFC_project.include_hemi_sorted,rois),ismember(app.dFC_project.include_hemi_sorted,rois));
    threshold = app.dFC_project.app_threshold;
    btwn_net_FC = 0;
    in_net_FC = 0;

    check = [167,201,209,217];
    num_net = sum(app.rois(check)>0);
    if num_net>1
        S.atlas = 1;
        nets = cell(4,1);
        nets{1} = get_network(S,1).include;
        nets{2} = get_network(S,3).include;
        nets{3} = get_network(S,2).include;
        nets{4} = get_network(S,4).include;

        color = zeros(numel(rois),3);
        colors = [0 0.4470 0.7410;
                  0.900 0.150 0.0980;
                  0.9290 0.6940 0.1250;
                  0.4940 0.1840 0.5560];
        count = 1;
        which_net = zeros(num_net,1);
        communities = zeros(numel(rois),1);
        for i=1:numel(check)
            if app.rois(check(i))>0
                color(ismember(rois,nets{i}),:) = repmat(colors(count,:),[sum(ismember(rois,nets{i})) 1]);
                which_net(count) = i;
                communities(ismember(rois,nets{i})) = count;
                count = count+1;
            end
        end
        try
            if num_net==2
                mygraph_2 = mygraph(ismember(rois,nets{which_net(1)}),ismember(rois,nets{which_net(2)}));
                btwn_net_FC = round(mean(mygraph_2,'all'),3);
            else
                btwn_net_FC = 0;
            end
        catch
            errordlg(lasterr)
        end
    elseif num_net==1 || num_net==0
        %     for i=1:max(comm)
        %         color(comm==i,:) = repmat(colors(i*floor(256/app.ncomm),:),[sum(comm==i) 1]);
        %     end
        % catch
            in_net_FC = round(mean(mygraph,'all'),3);
            colors = eval(app.plot_color); % problem with this not being defined?
            color = zeros(size(mygraph,1),3);
            A = mygraph;
        % try
        %     in_net_FC = round(mean(mygraph,'all'),3);
        %     colors = hsv;
        %     color = zeros(size(mygraph,1),3);
        %     comm = app.comm{state};
        %     ncomm = app.ncomm;
            A(A<0) = 0;
            rng_state = rng();
            rng(123)
            [comm,~] = community_louvain(A);
            if state~=1 && num_net==1
                A2 = app.dFC_project.state_dFC_mean{1}(ismember(app.dFC_project.include_hemi_sorted,rois),ismember(app.dFC_project.include_hemi_sorted,rois));
                A2(A2<0) = 0;
                [comm2,~] = community_louvain(A2);
                comm = harmonize_comm(comm2,comm);
            end
            rng(rng_state)
            for i=1:max(comm)
                if strcmp(app.plot_color,'default')
                    color(comm==i,:) = repmat(colors(i,:),[sum(comm==i) 1]);
                else
                    color(comm==i,:) = repmat(colors(i*floor(256/max(comm))-128/max(comm),:),[sum(comm==i) 1]);
                end
            end
        % end
    else
        color = 'intensity';
    end
elseif app.dFC_project.atlas==1 && nroi==86
    % app.rois
    % rois = 1:166;
    rois = app.dFC_project.include_hemi_sorted(ismember(app.dFC_project.include_hemi_sorted,rois));
    % if numel(rois)==166 && app.dFC_project.which_hemi==2
    %     rois = rois(1:78);
    %     hemi = 1:78;
    % elseif numel(rois)==166 && app.dFC_project.which_hemi==3
    %     rois = rois(89:166);
    %     hemi = 89:166;
    % % elseif numel(rois)==166
    % %     hemi = 1:166;
    % elseif app.dFC_project.which_hemi==2
    %     rois = rois(1:(numel(rois)/2));
    %     hemi = 1:78;
    % elseif app.dFC_project.which_hemi==3
    %     rois = rois((numel(rois)/2) + (1:(numel(rois)/2)));
    %     hemi = 89:166;
    % else
    %     hemi = 1:86;
    % end

    mygraph = app.dFC_project.state_dFC_mean{state};
    mygraph = mygraph(ismember(app.dFC_project.include_hemi_sorted,rois),ismember(app.dFC_project.include_hemi_sorted,rois));
    threshold = app.dFC_project.app_threshold;
    btwn_net_FC = 0;
    in_net_FC = 0;

    check = [167,201,209,217];
    num_net = sum(app.rois(check)>0);
    if num_net>1
        S.atlas = 1;
        nets = cell(4,1);
        nets{1} = get_network(S,1).include;
        nets{2} = get_network(S,3).include;
        nets{3} = get_network(S,2).include;
        nets{4} = get_network(S,4).include;

        color = zeros(numel(rois),3);
        colors = [0 0.4470 0.7410;
                  0.900 0.150 0.0980;
                  0.9290 0.6940 0.1250;
                  0.4940 0.1840 0.5560];
        count = 1;
        which_net = zeros(num_net,1);
        communities = zeros(numel(rois),1);
        for i=1:numel(check)
            if app.rois(check(i))>0
                color(ismember(rois,nets{i}),:) = repmat(colors(count,:),[sum(ismember(rois,nets{i})) 1]);
                which_net(count) = i;
                communities(ismember(rois,nets{i})) = count;
                count = count+1;
            end
        end
        try
            if num_net==2
                mygraph_2 = mygraph(ismember(rois,nets{which_net(1)}),ismember(rois,nets{which_net(2)}));
                btwn_net_FC = round(mean(mygraph_2,'all'),3);
            else
                btwn_net_FC = 0;
            end
        catch
            errordlg(lasterr)
        end
    elseif num_net==1 || num_net==0
        %     for i=1:max(comm)
        %         color(comm==i,:) = repmat(colors(i*floor(256/app.ncomm),:),[sum(comm==i) 1]);
        %     end
        % catch
            in_net_FC = round(mean(mygraph,'all'),3);
            colors = eval(app.plot_color);
            color = zeros(size(mygraph,1),3);
            A = mygraph;
        % try
        %     in_net_FC = round(mean(mygraph,'all'),3);
        %     colors = hsv;
        %     color = zeros(size(mygraph,1),3);
        %     comm = app.comm{state};
        %     ncomm = app.ncomm;
            % A(A<0) = 0;
            A = A - min(A,[],'all');
            rng_state = rng();
            rng(123)
            [comm,~] = community_louvain(A);
            if state~=1 && num_net==1
                A2 = app.dFC_project.state_dFC_mean{1}(ismember(app.dFC_project.include_hemi_sorted,rois),ismember(app.dFC_project.include_hemi_sorted,rois));
                % A2(A2<0) = 0;
                A2 = A2 - min(A2,[],'all');
                [comm2,~] = community_louvain(A2);
                comm = harmonize_comm(comm2,comm);
            end
            rng(rng_state)

            if num_net>0
                count = 1;
                which_net = zeros(num_net,1);
                for i=1:numel(check)
                    if app.rois(check(i))>0
                        which_net(count) = i;
                        count = count+1;
                    end
                end
                switch which_net
                    case 1
                        network = 'LN';
                    case 2
                        network = 'SN';
                    case 3
                        network = 'DAN';
                    case 4
                        network = 'DMN';
                    otherwise
                        network = '';
                end
                app.dFC_project.([network,'_communities']) = cell(max(comm),1);
            end

            for i=1:max(comm)
                if strcmp(app.plot_color,'default')
                    color(comm==i,:) = repmat(colors(i,:),[sum(comm==i) 1]);
                else
                    color(comm==i,:) = repmat(colors(i*floor(256/max(comm)-128/max(comm)),:),[sum(comm==i) 1]);
                end
                if num_net>0
                    app.dFC_project.([network,'_communities_state_',num2str(state)]){i} = app.dFC_project.all_roi_names(rois(comm==i));
                    disp([network,' state ',num2str(state),' community ',num2str(i),':'])
                    app.dFC_project.([network,'_communities_state_',num2str(state)]){i}
                end
            end
        % end
    else
        color = 'intensity';
    end
else
    mygraph = app.dFC_project.state_dFC_mean{state};
    if state~=1
        mygraph2 = app.dFC_project.state_dFC_mean{1};
    end

    in_net_FC = round(mean(mygraph,'all'),3);
    colors = eval(app.plot_color);
    color = zeros(size(mygraph,1),3);
    A = mygraph;
    A(A<0) = 0;
    rng_state = rng();
    rng(123)
    [comm,~] = community_louvain(A);
    rng(rng_state)
    
    ncomm = zeros(max(app.dFC_project.clusters),1);
    for i = 1:max(app.dFC_project.clusters)
        mygraph_temp = app.dFC_project.state_dFC_mean{i};
        A_temp = mygraph_temp;
        A_temp(A_temp<0) = 0;
        rng_state = rng();
        rng(123)
        [comm_temp,~] = community_louvain(A_temp);
        rng(rng_state)
        ncomm(i) = max(comm_temp);
    end

    max_ncomm = max(ncomm);

    if state~=1
        B = mygraph2;
        B(B<0) = 0;
        rng_state = rng();
        rng(123)
        [comm_B,~] = community_louvain(B);
        rng(rng_state)
        comm = harmonize_comm(comm_B,comm);
        for i=1:max(comm)
            if strcmp(app.plot_color,'default')
                color(comm==i,:) = repmat(colors(i,:),[sum(comm==i) 1]);
            else
                color(comm==i,:) = repmat(colors(i*floor(256/max_ncomm)-128/max_ncomm,:),[sum(comm==i) 1]);
            end
            app.dFC_project.(['communities_state_',num2str(state)]){i} = app.dFC_project.roi_names(comm==i);
            disp(['state ',num2str(state),' community ',num2str(i),':'])
            app.dFC_project.(['communities_state_',num2str(state)]){i}
        end
    else
        for i=1:max(comm)
            if strcmp(app.plot_color,'default')
                color(comm==i,:) = repmat(colors(i,:),[sum(comm==i) 1]);
            else
                color(comm==i,:) = repmat(colors(i*floor(256/max_ncomm)-128/max_ncomm,:),[sum(comm==i) 1]);
            end
            app.dFC_project.(['communities_state_',num2str(state)]){i} = app.dFC_project.roi_names(comm==i);
            disp(['state ',num2str(state),' community ',num2str(i),':'])
            app.dFC_project.(['communities_state_',num2str(state)]){i}
        end
    end

    threshold = app.dFC_project.app_threshold;
    rois = app.dFC_project.include_hemi_sorted;
    if app.dFC_project.which_hemi==2
        rois = rois(1:nroi);
    elseif app.dFC_project.which_hemi==3
        rois = rois(nroi + (1:nroi));
    end
    % color = 'intensity';
end



% 
% mygraph = app.dFC_project.state_dFC_mean{state};
% mygraph = mygraph(rois,rois);
% threshold = app.dFC_project.app_threshold;
% btwn_net_FC = 0;
% in_net_FC = 0;
% 
% check = [167,201,209,217];
% num_net = sum(app.rois(check)>0);
% if num_net>1
%     S.atlas = 1;
%     nets = cell(4,1);
%     nets{1} = get_network(S,1).include;
%     nets{2} = get_network(S,3).include;
%     nets{3} = get_network(S,2).include;
%     nets{4} = get_network(S,4).include;
%     color = zeros(numel(rois),3);
%     colors = [0 0.3 1;
%               1 0.3 0;
%               0.3 1 0;
%               0 1 0.3];
%     count = 1;
%     which_net = zeros(num_net,1);
%     for i=1:numel(check)
%         if app.rois(check(i))>0
%             color(ismember(rois,nets{i}),:) = repmat(colors(count,:),[sum(ismember(rois,nets{i})) 1]);
%             which_net(count) = i;
%             count = count+1;
%         end
%     end
%     try
%         if num_net==2
%             mean(mygraph,'all')
%             mygraph_2 = mygraph(ismember(rois,nets{which_net(1)}),ismember(rois,nets{which_net(2)}));
%             btwn_net_FC = round(mean(mygraph,'all'),3);
%         else
%             btwn_net_FC = 0;
%         end
%     catch
%         errordlg(lasterr)
%     end
% elseif num_net==1
%     % try
%     %     in_net_FC = round(mean(mygraph,'all'),3);
%     %     colors = hsv;
%     %     color = zeros(size(mygraph,1),3);
%     %     comm = app.comm{state};
%     %     ncomm = app.ncomm;
%     %     for i=1:max(comm)
%     %         color(comm==i,:) = repmat(colors(i*floor(256/app.ncomm),:),[sum(comm==i) 1]);
%     %     end
%     % catch
%         in_net_FC = round(mean(mygraph,'all'),3);
%         colors = hsv;
%         color = zeros(size(mygraph,1),3);
%         A = mygraph;
%         A(A<0) = 0;
%         rng_state = rng();
%         rng(123)
%         [comm,~] = community_louvain(A);
%         rng(rng_state)
%         for i=1:max(comm)
%             color(comm==i,:) = repmat(colors(i*floor(256/max(comm)),:),[sum(comm==i) 1]);
%         end
%     % end
% else
%     color = 'intensity';
% end
% 

atlas = 1;

if ~exist('atlas','var')
    atlas = 1;
end

t = rng;
rng(123)

mygraph_neg = mygraph;
mygraph(mygraph<0) = 0;
nroi = size(mygraph,1);
A = binarize(mygraph,threshold);

if atlas==2
    vox_coor = get_schaefer_neurosynth_vox_coor();
else
    vox_coor = get_AAL3_vox_coor();
end

vox_coor = vox_coor(rois,:);

% size(vox_coor)
% numel(rois)

% S = struct(brush=0, dark=false, detail=0, grid = false);
% spm_glass(ones(numel(rois),1),vox_coor,S);
hold(app.UIAxes,'on')

if ischar(color) && strcmp(color,'intensity')
    colors = 1;
    % unique_colors = turbo;
    unique_colors = zeros(256,3);
    unique_colors(1:128,2) = flip([1:128]/128);
    unique_colors(1:128,1) = flip([1:128]/128);
    unique_colors(1:128,3) = 1;
    unique_colors(129:256,3) = flip(0.5+0.5*[1:128]/128);
elseif ischar(color) && strcmp(color,'intensity_neg')
    colors = 2;
    unique_colors = zeros(256,3);
    % unique_colors(65:128,2) = [1:64]/64;
    % unique_colors(129:192) = flip([1:64]/64);
    % unique_colors(1:64,3) = 0.5+0.5*[1:64]/64;
    % unique_colors(65:192,3) = flip([1:128]/128);
    % unique_colors(65:192,1) = [1:128]/128;
    % unique_colors(193:256,1) = flip(0.5+0.5*[1:64]/64);
    unique_colors(129:256,1) = 1;
    unique_colors(129:256,2) = (127-(0:127))/127;
    unique_colors(129:256,3) = (127-(0:127))/127;
    unique_colors(1:128,1) = (0:127)/127;
    unique_colors(1:128,2) = (0:127)/127;
    unique_colors(1:128,3) = 1;
elseif ischar(color)
    disp("'color' argument must be 'intensity' or a matrix of RGB triplets")
else
    colors = 0;
    unique_colors = unique(color,'rows');
end

max_weight = max(mygraph,[],2);
for i=1:nroi
    if max(mygraph(i,[1:nroi]~=i))<threshold
        max_weight(i) = -1;
    end
end

test = zeros(length(max_weight),size(unique_colors,1));

% if colors~=2
%     for ii=1:size(unique_colors,1)
%         if colors==1
%             ii_color = floor(256*(max_weight-threshold)/(max(max_weight))-threshold)==ii;
%             face_alpha = log(ii)/log(256);
%         else
%             ii_color = sum(color==unique_colors(ii,:),2)==3;
%             face_alpha = 1;
%         end
%         xy_axial = mni2xy_v2(vox_coor(ii_color,1),vox_coor(ii_color,2),vox_coor(ii_color,3),'axial') + [-80,95];
%         xy_coronal = mni2xy_v2(vox_coor(ii_color,1),vox_coor(ii_color,2),vox_coor(ii_color,3),'coronal');
%         xy_sagittal = mni2xy_v2(vox_coor(ii_color,1),vox_coor(ii_color,2),vox_coor(ii_color,3),'sagittal') + [-12,0];
%         scatter(app.UIAxes,xy_axial(:,1),xy_axial(:,2),12,...
%             'MarkerFaceColor',unique_colors(ii,:),...
%             'MarkerEdgeColor','none',...
%             'MarkerFaceAlpha',face_alpha)
%         scatter(app.UIAxes,xy_coronal(:,1),xy_coronal(:,2),12,...
%             'MarkerFaceColor',unique_colors(ii,:),...
%             'MarkerEdgeColor','none',...
%             'MarkerFaceAlpha',face_alpha)
%         scatter(app.UIAxes,xy_sagittal(:,1),xy_sagittal(:,2),12,...
%             'MarkerFaceColor',unique_colors(ii,:),...
%             'MarkerEdgeColor','none',...
%             'MarkerFaceAlpha',face_alpha)
%     end
% else
%     for i=1:nroi
%         xy_axial = mni2xy_v2(vox_coor(i,1),vox_coor(i,2),vox_coor(i,3),'axial') + [-80,95];
%         xy_coronal = mni2xy_v2(vox_coor(i,1),vox_coor(i,2),vox_coor(i,3),'coronal');
%         xy_sagittal = mni2xy_v2(vox_coor(i,1),vox_coor(i,2),vox_coor(i,3),'sagittal') + [-12,0];
%         scatter(app.UIAxes,xy_axial(:,1),xy_axial(:,2),8,...
%             'MarkerFaceColor',[0.5 0.5 0.5],...
%             'MarkerEdgeColor','none')
%         scatter(app.UIAxes,xy_coronal(:,1),xy_coronal(:,2),8,...
%             'MarkerFaceColor',[0.5 0.5 0.5],...
%             'MarkerEdgeColor','none')
%         scatter(app.UIAxes,xy_sagittal(:,1),xy_sagittal(:,2),8,...
%             'MarkerFaceColor',[0.5 0.5 0.5],...
%             'MarkerEdgeColor','none')
%     end
% end

% stop= 0;
for i = 1:size(vox_coor,1)
    for j=1:size(vox_coor,1)
        % if and(node_coor(i,1)==0, node_coor(i,2)==0) && stop==0
        %     disp(['i = ', num2str(i)])
        %     disp(['j = ', num2str(j)])
        %     stop = 1;
        % end
        if (i~=j && A(i,j)==1) || (i~=j && colors==2 && abs(mygraph_neg(i,j))>threshold)
            xy_axial = mni2xy_v2([vox_coor(i,1),vox_coor(j,1)],[vox_coor(i,2),vox_coor(j,2)],[vox_coor(i,3),vox_coor(j,3)],'axial');
            xy_sagittal = mni2xy_v2([vox_coor(i,1),vox_coor(j,1)],[vox_coor(i,2),vox_coor(j,2)],[vox_coor(i,3),vox_coor(j,3)],'sagittal');
            xy_coronal = mni2xy_v2([vox_coor(i,1),vox_coor(j,1)],[vox_coor(i,2),vox_coor(j,2)],[vox_coor(i,3),vox_coor(j,3)],'coronal');
            
            xy_axial = xy_axial + [-80,95];
            xy_sagittal = xy_sagittal + [-12,0];

            if ~ischar(color) && size(color,1)>1
                patchline(xy_axial(:,1),xy_axial(:,2),...
                    app.UIAxes,...
                    'linewidth',1.5,...
                    'edgealpha',(mygraph(i,j)-threshold)/(max(mygraph,[],'all')-threshold),...
                    'edgecolor',mean([color(j,:);color(i,:)],1))
    
                patchline(xy_sagittal(:,1),xy_sagittal(:,2),...
                    app.UIAxes,...
                    'linewidth',1.5,...
                    'edgealpha',(mygraph(i,j)-threshold)/(max(mygraph,[],'all')-threshold),...
                    'edgecolor',mean([color(j,:);color(i,:)],1))
    
                patchline(xy_coronal(:,1),xy_coronal(:,2),...
                    app.UIAxes,...
                    'linewidth',1.5,...
                    'edgealpha',(mygraph(i,j)-threshold)/(max(mygraph,[],'all')-threshold),...
                    'edgecolor',mean([color(j,:);color(i,:)],1))
            elseif ~ischar(color)
                patchline(xy_axial(:,1),xy_axial(:,2),...
                    app.UIAxes,...
                    'linewidth',1.5,...
                    'edgealpha',(mygraph(i,j)-threshold)/(max(mygraph,[],'all')-threshold),...
                    'edgecolor',color)
    
                patchline(xy_sagittal(:,1),xy_sagittal(:,2),...
                    app.UIAxes,...
                    'linewidth',1.5,...
                    'edgealpha',(mygraph(i,j)-threshold)/(max(mygraph,[],'all')-threshold),...
                    'edgecolor',color)
    
                patchline(xy_coronal(:,1),xy_coronal(:,2),...
                    app.UIAxes,...
                    'linewidth',1.5,...
                    'edgealpha',(mygraph(i,j)-threshold)/(max(mygraph,[],'all')-threshold),...
                    'edgecolor',color)
            elseif ~strcmp(color,'intensity') && ~strcmp(color,'intensity_neg')
                disp("'color' argument must be 'intensity' or a matrix of RGB triplets")
                return
            elseif strcmp(color,'intensity_neg')
                edge_color = 128 + floor(32*mygraph_neg(i,j));
                % disp(edge_color)
                edge_color(edge_color>256) = 256;
                patchline(xy_axial(:,1),xy_axial(:,2),...
                    app.UIAxes,...
                    'linewidth',1.5,...
                    'edgealpha',(abs(mygraph_neg(i,j))-threshold)/(max(abs(mygraph_neg),[],'all')-threshold),...
                    'edgecolor',unique_colors(edge_color,:))
    
                patchline(xy_sagittal(:,1),xy_sagittal(:,2),...
                    app.UIAxes,...
                    'linewidth',1.5,...
                    'edgealpha',(abs(mygraph_neg(i,j))-threshold)/(max(abs(mygraph_neg),[],'all')-threshold),...
                    'edgecolor',unique_colors(edge_color,:))
    
                patchline(xy_coronal(:,1),xy_coronal(:,2),...
                    app.UIAxes,...
                    'linewidth',1.5,...
                    'edgealpha',(abs(mygraph_neg(i,j))-threshold)/(max(abs(mygraph_neg),[],'all')-threshold),...
                    'edgecolor',unique_colors(edge_color,:))
            else
                edge_color = max(floor((mygraph(i,j)-threshold)*(256/(1.5-threshold))),1);
                edge_color(edge_color>256) = 256;
                patchline(xy_axial(:,1),xy_axial(:,2),...
                    app.UIAxes,...
                    'linewidth',1.5,...
                    'edgealpha',(mygraph(i,j)-threshold)/(max(mygraph,[],'all')-threshold),...
                    'edgecolor',unique_colors(edge_color,:))
    
                patchline(xy_sagittal(:,1),xy_sagittal(:,2),...
                    app.UIAxes,...
                    'linewidth',1.5,...
                    'edgealpha',(mygraph(i,j)-threshold)/(max(mygraph,[],'all')-threshold),...
                    'edgecolor',unique_colors(edge_color,:))
    
                patchline(xy_coronal(:,1),xy_coronal(:,2),...
                    app.UIAxes,...
                    'linewidth',1.5,...
                    'edgealpha',(mygraph(i,j)-threshold)/(max(mygraph,[],'all')-threshold),...
                    'edgecolor',unique_colors(edge_color,:))
            end
        end
    end
end

% pbaspect([3 1 1])

I = imread('brain_background_v2_inline.png');
% colors = gray;
% colormap(flip(colors,1))
colormap(app.UIAxes,'gray')
I = double(I);
I = I + 256-max(I,[],'all');
h = image(app.UIAxes,[-80 184] ,[200 100],I); 
uistack(h,'bottom')

networks = {
    'Language'
    'Dorsal Attention'
    'Salience'
    'Default Mode'};

colors = [0 0.4470 0.7410;
          0.900 0.150 0.0980;
          0.9290 0.6940 0.1250;
          0.4940 0.1840 0.5560];

if exist('num_net','var') && num_net==2
    text(app.UIAxes,100,195,['Between-network FC = ',num2str(btwn_net_FC)])
    text(app.UIAxes,150,100+5*num_net,networks{which_net(1)},'Color',colors(1,:),'FontSize',9)
    text(app.UIAxes,150,100+5*(num_net-1),networks{which_net(2)},'Color',colors(2,:),'FontSize',9)
elseif exist('num_net','var') && num_net==1
    text(app.UIAxes,100,195,['Within-network FC = ',num2str(in_net_FC)])
elseif exist('num_net','var')
    for i=1:num_net
        text(app.UIAxes,150,100+5*(num_net+1-i),networks{which_net(i)},'Color',colors(i,:),'FontSize',9)
    end
else
    text(app.UIAxes,120,195,['Within-network FC = ',num2str(in_net_FC)])
end
hold(app.UIAxes,'off')

if numel(rois)>0
    try
        plot_BU_graph_comm_app
    catch
        errordlg(lasterr)
    end
end

% xlim([-80 180])
% ylim([100 200])
app.UIAxes.XLim = [-80 180];
app.UIAxes.YLim = [100 200];

% switch orient
%     case 'axial'
%         xlim([0 115])
%         ylim([95 180])
%         set(gca,'ydir','reverse')
%         view([90 -90])
%     case 'sagittal'
%         xlim([0 105])
%         ylim([0 95])
%     case 'coronal'
%         xlim([105 210])
%         ylim([5 95])
%         set(gca,'xdir','reverse')
%     otherwise
% end

rng(t)