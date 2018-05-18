function ee = DatPFC_PhotonCounting(ee,XML_Settings,Debug)
%% Excerp from PulseQuantities_PhotonCounting
if ~exist('Debug','var')
    Debug = 0;
end
%% Initialize variables

pmt_chs = 1:122;

per_channel_event_size = [XML_Settings.max_num_pulses length(pmt_chs) 1];

ee.rqs.peak_width_samples                   = nan(per_channel_event_size);
ee.rqs.spike_count                          = nan(per_channel_event_size);
% new RQs --pfs
%% Compute RQs

try
    loop_max = min([sum(isfinite(ee.rqs.index_kept_sumpods(:,1))) XML_Settings.max_num_pulses]);
    if loop_max == 0
        if Debug
            disp('Fail to run photon counting.')
        end
        ee.info.PhotonCountSuccess = 0;
        ee.info.PhotonCountError = 'No pulse found.';
        return
    end
    % for every pulse (NOT pod)
    for pp = 1:loop_max
        for ch = pmt_chs
            % Make slice of data
            peak_cut = inrange(ee.ch(ch).pod_time_samples,ee.rqs.pulse_start_samples(pp,1)-1,ee.rqs.pulse_end_samples(pp,1)+1);
            if ~isempty(peak_cut) && sum(peak_cut) > 0
                if strcmp(XML_Settings.height_units, 'mV') == 1                          
                    peak_data_height = ee.ch(ch).pod_data_mV(peak_cut);
                elseif strcmp(XML_Settings.height_units, 'phe') == 1
                    peak_data_height = ee.ch(ch).pod_data_phe_per_sample(peak_cut);
                end
                peak_time_samples = ee.ch(ch).pod_time_samples(peak_cut);

                % This is to find which baseline_mV to use. Since
                % we are cutting pods in time, it's not
                % straightforward to know which pod baseline 'twas                            
                % But this statement takes care of that

                %Commented out since unused
                %peak_inds = unique(cumsum(ismember(ee.ch(ch).pod_time_samples,ee.ch(ch).pod_start_samples)));

                % new RQs related with the photon couting module
                %test 0
                peak_above_threshold = peak_time_samples(peak_data_height>XML_Settings.threshold);

                if ~isempty(peak_above_threshold)
                    if numel(peak_data_height) > 2
                        ee.rqs.peak_width_samples(pp,ch,1) = max(peak_above_threshold)-min(peak_above_threshold)+1;
                        peak_height_phe = peak_data_height > XML_Settings.threshold;
                        Variacao = (peak_data_height(2:end)-peak_data_height(1:end-1)).*peak_height_phe(1:end-1);
                        ee.rqs.spike_count(pp,ch,1) = sum(Variacao(2:end)<0 & Variacao(1:end-1)>=0);
                    end
                end

            end % fi
        end % for ch

    end % for pulse pp
    ee.info.PhotonCountSuccess = 1;
    ee.info.PhotonCountError = '';
    
catch exception
    if Debug
        disp('Fail to run photon counting.')
    end
    ee.info.PhotonCountSuccess = 0;
    ee.info.PhotonCountError = exception.identifier;
    return
end
