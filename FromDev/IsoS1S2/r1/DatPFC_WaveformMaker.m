function ee = DatPFC_WaveformMaker(raw,TrgTime,EventWindow,SR_Version)
    try
        if strcmp(SR_Version,'SR2.1')
            TimeOffset = 0;
        elseif strcmp(SR_Version,'SR2.0')
            TimeOffset = -24;
        else
            TimeOffset = -24;
        end
        if ~isfield(raw,'PodTimeList')
            raw.PodTimeList = cell(128,1);
            for ch = 1:min(128,length(raw.ch))
                if ~isempty(raw.ch(ch).pod)
                    raw.PodTimeList{ch,1} = vertcat(raw.ch(ch).pod.timestamp);
                end
            end
        end
        %% make sumpod
        ee.ch_map = 1:122;
        ee.info.EventWindow = EventWindow;
        ee.timestamp = TrgTime;
        ee.empty = 0;
        ee.thr = 0;
        ee.trigger_pulse_timestamp_samples = [];
        ee.ch = struct(...
            'empty',cell(1,122),...
            'pod_start_samples',cell(1,122),...
            'pod_length_samples',cell(1,122),...
            'pod_baseline',cell(1,122),...
            'pod_baseline_mV',cell(1,122),...
            'pod_data',cell(1,122),...
            'pod_time_samples',cell(1,122),...
            'pod_data_mV',cell(1,122)...
            );

        PodInEvtList = find((EventWindow(1)+TrgTime <= raw.PodTimeList{128,1}) & (raw.PodTimeList{128,1} <= EventWindow(2)+TrgTime));
        if ~isempty(PodInEvtList)
            ee.trigger_pulse_timestamp_samples = vertcat(raw.ch(128).pod.timestamp);
        end
        for ch = 1:122
            PodInEvtList = find((EventWindow(1)+TrgTime <= raw.PodTimeList{ch,1}) & (raw.PodTimeList{ch,1} <= EventWindow(2)+TrgTime));
            if ~isempty(PodInEvtList)
                ee.ch(ch).empty = 0;
                ee.ch(ch).pod_start_samples = vertcat(raw.ch(ch).pod(PodInEvtList).timestamp ) - TrgTime + TimeOffset; 
                ee.ch(ch).pod_length_samples = vertcat(raw.ch(ch).pod(PodInEvtList).length);
                ee.ch(ch).pod_data = vertcat(raw.ch(ch).pod(PodInEvtList).pod_data);

                ee.ch(ch).pod_baseline = zeros(length(PodInEvtList),1);
                ee.ch(ch).pod_baseline_mV = zeros(length(PodInEvtList),1);
                ee.ch(ch).pod_data_mV = [];
                ee.ch(ch).pod_time_samples = [];
                for ps = 1:length(PodInEvtList)
                    ee.ch(ch).pod_baseline(ps,1) = raw.ch(ch).pod(PodInEvtList(ps)).baseline;
                    ee.ch(ch).pod_baseline_mV(ps,1) = -(mean(double(raw.ch(ch).pod(PodInEvtList(ps)).pod_data(1:22)))*2000/(2^14) - 1900) ;
                    ee.ch(ch).pod_time_samples = horzcat(ee.ch(ch).pod_time_samples,(raw.ch(ch).pod(PodInEvtList(ps)).timestamp - TrgTime + (0:raw.ch(ch).pod(PodInEvtList(ps)).length-1)) + TimeOffset);
                    ee.ch(ch).pod_data_mV = vertcat(ee.ch(ch).pod_data_mV,-(raw.ch(ch).pod(PodInEvtList(ps)).pod_data * 2000/(2^14) - 1900)- ee.ch(ch).pod_baseline_mV(ps,1));
                end
            else
                ee.ch(ch).empty = 1;

            end
        end
        ee.info.EventConstructionSuccess = 1;
        ee.info.EventConstructionError = '';
    catch exception
        ee.info.EventConstructionSuccess = 0;
        ee.info.EventConstructionError = exception.identifier;
    end