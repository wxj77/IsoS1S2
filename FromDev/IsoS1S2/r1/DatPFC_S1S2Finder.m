function ee = DatPFC_S1S2Finder(ee,XML_Settings,Debug)
if ~exist('Debug','var')
    Debug = 0;
end



%% Excerp/Adapted from S1S2Pairing_Naive.m

pulse_event_size = [XML_Settings.max_num_pulses 1];


MaxDriftLength = 50000; %In samples

ee.rqs.NumberOfLeadingS1 = nan(pulse_event_size);
ee.rqs.NumberOfLeadingS2 = nan(pulse_event_size);
ee.rqs.NumberOfLeadingS3 = nan(pulse_event_size);
ee.rqs.NumberOfLeadingS4 = nan(pulse_event_size);
ee.rqs.NumberOfLeadingS5 = nan(pulse_event_size);

ee.rqs.NumberOfTailingS1 = nan(pulse_event_size);
ee.rqs.NumberOfTailingS2 = nan(pulse_event_size);
ee.rqs.NumberOfTailingS3 = nan(pulse_event_size);
ee.rqs.NumberOfTailingS4 = nan(pulse_event_size);
ee.rqs.NumberOfTailingS5 = nan(pulse_event_size);

ee.rqs.s1s2_pairing      = nan(pulse_event_size);
ee.rqs.z_drift_samples   = nan(pulse_event_size);

if XML_Settings.max_num_pulses == 0
    if Debug
        disp('Fail to run S1-S2 pairing.')
    end
    ee.info.S1S2PairingSuccess = 0;
    ee.info.S1S2PairingError = 'No pulse found.';
    return
end

if (sum(ee.rqs.pulse_classification == 1) < 1)
    ee.info.S1S2PairingSuccess = 0;
    ee.info.S1S2PairingError = 'This dat file has no S1. Can not progress to pair the pulses.';
    return
end
if (sum(ismember(ee.rqs.pulse_classification,[2 4])) < 1)
    ee.info.S1S2PairingSuccess = 0;
    ee.info.S1S2PairingError = 'This dat file has no S2. Can not progress to pair the pulses.';
    return
end


%% Loop per event and assign S1 S2 pairs
try
    for ii_pp = 1:XML_Settings.max_num_pulses
        for ii_pulsetype = 1:5
            ee.rqs.(['NumberOfLeadingS' num2str(ii_pulsetype)])(ii_pp) = ...
                sum(...
                ( (ee.rqs.aft_t0_samples - ee.rqs.aft_t0_samples(ii_pp)) < 0) & ...
                ( (ee.rqs.aft_t0_samples - ee.rqs.aft_t0_samples(ii_pp)) > -MaxDriftLength) & ...
                (ee.rqs.pulse_classification == ii_pulsetype)...
                );
            ee.rqs.(['NumberOfTailingS' num2str(ii_pulsetype)])(ii_pp) = ...
                sum(...
                ( (ee.rqs.aft_t0_samples - ee.rqs.aft_t0_samples(ii_pp)) > 0) & ...
                ( (ee.rqs.aft_t0_samples - ee.rqs.aft_t0_samples(ii_pp)) < MaxDriftLength) & ...
                (ee.rqs.pulse_classification == ii_pulsetype)...
                );
        end
    end

%% Idea

%1. Each S1 gets an ID from 1 to total number of S1
%2. Each S2/SE is paired with a specific S1 will get the same ID as the S1

%DP module always selects the first S1 in an event.
    
    % find the indices for S1 and S2 pulses
    s1_inds = find(ee.rqs.pulse_classification == 1);
    s2_inds = find( (ee.rqs.pulse_classification == 2) | (ee.rqs.pulse_classification == 4) );
    ee.rqs.s1s2_pairing(s1_inds) = 1:length(s1_inds);
    
    if ~isempty(s1_inds) && ~isempty(s2_inds)
        for ii = 1:length(s2_inds)
            pps2 = s2_inds(ii);
            if any((ee.rqs.pulse_classification == 1) & ( (ee.rqs.aft_t0_samples(pps2) - ee.rqs.aft_t0_samples) < MaxDriftLength ) & ( (ee.rqs.aft_t0_samples(pps2) - ee.rqs.aft_t0_samples) > 0 ))
                ee.rqs.s1s2_pairing(pps2) = ...
                    ee.rqs.s1s2_pairing(find( (ee.rqs.pulse_classification == 1) & ( (ee.rqs.aft_t0_samples(pps2) - ee.rqs.aft_t0_samples) < MaxDriftLength )& ( (ee.rqs.aft_t0_samples(pps2) - ee.rqs.aft_t0_samples) > 0 ),1,'first')); 
                ee.rqs.z_drift_samples(pps2) = ee.rqs.aft_t0_samples(pps2) - ...
                    ee.rqs.aft_t0_samples(find( (ee.rqs.pulse_classification == 1) & ( (ee.rqs.aft_t0_samples(pps2) - ee.rqs.aft_t0_samples) < MaxDriftLength )& ( (ee.rqs.aft_t0_samples(pps2) - ee.rqs.aft_t0_samples) > 0 ),1,'first')); 
            end
            %{
            ee.rqs.s1s2_pairing_ClosestS1(pps2) = ...
                ee.rqs.s1s2_pairing_ClosestS1(find( (ee.rqs.pulse_classification == 1) & ( (ee.rqs.aft_t0_samples(pps2) - ee.rqs.aft_t0_samples) < MaxDriftLength ),1,'last')); 
            ee.rqs.z_drift_samples_ClosestS1(pps2) = ee.rqs.aft_t0_samples(pps2) - ...
                ee.rqs.aft_t0_samples(find( (ee.rqs.pulse_classification == 1) & ( (ee.rqs.aft_t0_samples(pps2) - ee.rqs.aft_t0_samples) < MaxDriftLength ),1,'last')); 
            %}
        end
        
    end
    ee.info.S1S2PairingSuccess = 1;
    ee.info.S1S2PairingError = '';
catch exception
    if Debug
        disp('Fail to run S1-S2 pairing.')
    end
    ee.info.S1S2PairingSuccess = 0;
    ee.info.S1S2PairingError = exception.identifier;
end





