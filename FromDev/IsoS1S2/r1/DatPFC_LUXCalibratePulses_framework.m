function ee = DatPFC_LUXCalibratePulses_framework(ee, mVns_per_phe, amp_gain, Debug)
if ~exist('Debug','var')
    Debug = 0;
end
try
    if isfield(ee,'ch_map')
        %ch_map = event_struct(1).ch_map;
        ch_map = 1:122; % this should not be hard-coded, but the previous line is not good either!
        %fprintf('warning: ch_map is hard-coded (LUXCalibratePulses_framework.m)\n');
    else
        fprintf('channel map not defined in event_struct\n');
        fprintf('assuming channel map is 1:number of channels\n');
        ch_map = 1:length(ee.ch);
    end

    conversion_factor_mV_to_phe_per_10ns = 10 ./(mVns_per_phe .* amp_gain) ;
    if ~ee.empty
        for ii = 1:length(ch_map)
            ch = ch_map(ii);
            if ee.ch(ch).empty == 0 % if it's not empty
                data = ee.ch(ch).pod_data_mV;
                data(data<ee.thr & data>-ee.thr) = 0; % symmetric threshold always applied
                ee.ch(ch).pod_data_phe_per_sample = ...
                    data ...
                    .* conversion_factor_mV_to_phe_per_10ns(ii);
            end
        end
    end
    ee.info.PulseCalibrationSuccess = 1;
    ee.info.PulseCalibrationError = '';
catch exception
    if Debug
        disp('Fail to run pulse_calibration.')
    end
    for ch = 1:122
        ee.ch(ch).pod_data_phe_per_sample = [];
    end
    ee.info.PulseCalibrationSuccess = 0;
    ee.info.PulseCalibrationError = exception.identifier;
end


