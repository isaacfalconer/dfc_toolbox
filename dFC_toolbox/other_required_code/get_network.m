function FC_project = get_network(FC_project,network)

old_path = pwd;

switch FC_project.atlas
    case 1
        LN = [1,2,3,4,5,6,7,8,9,10,11,12,15,16,25,26,57,58,61,62,65,66,67,68,81,82,83,84,85,86,87,88,89,90];
        SN = [33,34,147:152];
        DAN = [1,2,5,6,51,52,61,62];
        DMN = [19,20,25,26,37,38,63,64,69,70];
        LN_alt = [3,4,7,8,9,10,11,12,15,16,57,58,65,66,67,68,81,82,83,84,85,86,87,88,89,90];
        WC =    [1	2	3	4	5	6	7	8	9	10	11	12	13	14	15	16 ...
        	    17	18	19	20	21	22	23	24	25	26	27	28	29	30	31	32 ...
                33	34	35	36	37	38	41	42	45	46	47	48	49	50	51	52 ...
                53	54	55	56	57	58	59	60	61	62	63	64	65	66	67	68 ...
                69  70	71	72	79	80	81	82	83	84	85	86	87	88	89	90 ...
                147 148	149	150	151	152];

        switch network
            case 1
                FC_project.include = LN;
            case 2
                FC_project.include = SN;
            case 3
                FC_project.include = DAN;
            case 4
                FC_project.include = DMN;
            case 5
                nroi = str2double(inputdlg('How many ROIs would you like to include?'));
                FC_project.include = 2*randsample(166/2,nroi/2);
                FC_project.include = [FC_project.include;FC_project.include-1];
                FC_project.include = sort(FC_project.include);
            case 7
                FC_project.include = WC;
            otherwise
                disp oops!
        end
    case 2
        all_networks = get_schaefer_networks(FC_project);
        networks = all_networks(network);
        cd(FC_project.CONN_results_dir)
        names = schaefer_roi_names();
        k = 1;
        for i=1:numel(names)
            split_name = split(names{i},'_');
            net_name = split_name{3};
            for j=1:numel(networks)
                if strcmp(net_name,networks(j))
                    include(k) = i;
                    k = k+1;
                end
            end
        end
        FC_project.include = include;
    case 3
        all_networks = get_schaefer_networks_200(FC_project);
        networks = all_networks(network);
        cd(FC_project.CONN_results_dir)
        names = schaefer_roi_names_200();
        k = 1;
        for i=1:numel(names)
            split_name = split(names{i},'_');
            net_name = split_name{3};
            for j=1:numel(networks)
                if strcmp(net_name,networks(j))
                    include(k) = i;
                    k = k+1;
                end
            end
        end
        FC_project.include = include; 
    case 4
        LN = [50:78];
        SN = [79:86];
        DAN = [36:49];
        DMN = [1:35];

        switch network
            case 1
                FC_project.include = LN;
            case 2
                FC_project.include = SN;
            case 3
                FC_project.include = DAN;
            case 4
                FC_project.include = DMN;
            case 5
                nroi = str2double(inputdlg('How many ROIs would you like to include?'));
                FC_project.include = 2*randsample(166/2,nroi/2);
                FC_project.include = [FC_project.include;FC_project.include-1];
                FC_project.include = sort(FC_project.include);
            otherwise
                disp oops!
        end
    otherwise
        disp oops!
end

cd(old_path)

end