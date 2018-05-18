function ee = DatPFC_LUXSumPOD_framework(ee,ch_map,SR_Version,Debug)

if ~exist('SR_Version','var')
    disp('SR_Version is not specified. Use default SR2.0.')
    SR_Version = 'SR2.0';
end    
if ~exist('Debug','var')
    Debug = 0;
end
try
    if strcmp(SR_Version,'SR2.0')
        % Shortcut
        ch_map = 1:122;

        % Only get index of events that are not flagged as empty
        good_evts = find(cell2mat({ee.empty}) == 0);

        % For each event with something in it
        for evt = good_evts

            data = cat(1,ee(evt).ch(ch_map).pod_data_phe_per_sample);
            % DO NOT USE THIS (it would need to be implemented per channel) -pfs
            %data(data<ee(1).thr & data>-ee(1).thr) = 0; % symmetric threshold always applied
            times = cat(2,ee(evt).ch(ch_map).pod_time_samples)';

            inds = times - min(times) + 1;

            times_with_data = unique(times);
            inds_with_data = times_with_data - min(times) + 1;

            sumpod_data_phe_per_sample = accumarray(inds,data);
            sumpod_time_samples = min(times):max(times);

            ee(evt).sumpod_data_phe_per_sample = sumpod_data_phe_per_sample(inds_with_data);
            ee(evt).sumpod_time_samples = sumpod_time_samples(inds_with_data);

            %if 1 % additional code to make thresholded sumpod -pfs
            %	data_thr = data;
            %	data_thr(data_thr<ee(1).thr & data_thr>-ee(1).thr) = 0; % symmetric threshold
            %	sumpod_data_thr_phe_per_sample = accumarray(inds,data_thr);
            %	ee(evt).sumpod_data_thr_phe_per_sample = sumpod_data_thr_phe_per_sample(inds_with_data);
            %end

            % Sanity check
            if 0
                figure(432); clf; 
                h1 = subplot(2,1,1); 
                plot(ee(evt).sumpod_time_samples,ee(evt).sumpod_data_phe_per_sample,'k.-');

                h2 = subplot(2,1,2);
                plot(times,data,'.-');

                linkaxes([h1 h2],'xy')
                keyboard
            end

            ind_edges = find(diff(inds_with_data) > 1)';
            sumpod_end_samples =  ee(evt).sumpod_time_samples([ind_edges end]);

            ee(evt).sumpod_start_samples = ee(evt).sumpod_time_samples([1 ind_edges+1]);
            ee(evt).sumpod_length_samples = sumpod_end_samples - ee(evt).sumpod_start_samples + 1;

        end % evt
        ee(evt).info.SumPODSuccess = 1;
        ee(evt).info.SumPODError = '';
        ee(evt).info.SumPODVersion = SR_Version;
    elseif strcmp(SR_Version,'SR2.1')
        % Shortcut
        if ~exist('ch_map','var')
            ch_map = 1:122;
        end

        % Only get index of events that are not flagged as empty
        good_evts = find(cell2mat({ee.empty}) == 0);

        % For each event with something in it
        for evt = good_evts
           if mean([ee(evt).ch(ch_map).empty]) < 1
            if sum(~[ee(evt).ch(ch_map).empty]) > 0
            data = cat(1,ee(evt).ch(ch_map).pod_data_phe_per_sample);
            % DO NOT USE THIS (it would need to be implemented per channel) -pfs
            %data(data<ee(1).thr & data>-ee(1).thr) = 0; % symmetric threshold always applied
            times = cat(2,ee(evt).ch(ch_map).pod_time_samples)';

            inds = times - min(times) + 1;

            times_with_data = unique(times);
            inds_with_data = times_with_data - min(times) + 1;

            sumpod_data_phe_per_sample = accumarray(inds,data);
            sumpod_time_samples = min(times):max(times);

            ee(evt).sumpod_data_phe_per_sample = sumpod_data_phe_per_sample(inds_with_data);
            ee(evt).sumpod_time_samples = sumpod_time_samples(inds_with_data);

            %if 1 % additional code to make thresholded sumpod -pfs
            %	data_thr = data;
            %	data_thr(data_thr<ee(1).thr & data_thr>-ee(1).thr) = 0; % symmetric threshold
            %	sumpod_data_thr_phe_per_sample = accumarray(inds,data_thr);
            %	ee(evt).sumpod_data_thr_phe_per_sample = sumpod_data_thr_phe_per_sample(inds_with_data);
            %end

            % Sanity check
            if 0
                figure(432); clf; 
                h1 = subplot(2,1,1); 
                plot(ee(evt).sumpod_time_samples,ee(evt).sumpod_data_phe_per_sample,'k.-');

                h2 = subplot(2,1,2);
                plot(times,data,'.-');

                linkaxes([h1 h2],'xy')
                keyboard
            end

            ind_edges = find(diff(inds_with_data) > 1)';
            sumpod_end_samples =  ee(evt).sumpod_time_samples([ind_edges end]);

            ee(evt).sumpod_start_samples = ee(evt).sumpod_time_samples([1 ind_edges+1]);
            ee(evt).sumpod_length_samples = sumpod_end_samples - ee(evt).sumpod_start_samples + 1;
            end
           end
        end % evt
        ee(evt).info.SumPODSuccess = 1;
        ee(evt).info.SumPODError = '';
        ee(evt).info.SumPODVersion = SR_Version;
    else
        error('Do not recognize SR version.');
    end
catch exception
    if Debug
        disp('Fail to run LUXSumPOD.')
    end
    % Only get index of events that are not flagged as empty
    good_evts = find(cell2mat({ee.empty}) == 0);

    % For each event with something in it
    for evt = good_evts
        ee(evt).sumpod_data_phe_per_sample = [];
        ee(evt).sumpod_time_samples = [];
        ee(evt).sumpod_start_samples = [];
        ee(evt).sumpod_length_samples = [];

        ee(evt).info.SumPODSuccess = 0;
        ee(evt).info.SumPODError = exception.identifier;
    end
end