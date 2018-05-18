function timing = DatPFC_PulseTiming_HeightTiming_Converted_HeightTiming(pulse_data_phe,threshold)
    timing = nan(8,1);
    last = length(pulse_data_phe)-1;
    
    %Find peak value and location (t1)
    [peak_value, t1_idx_temp] = max(pulse_data_phe);
    t1_idx = t1_idx_temp -1;
    timing(4) = t1_idx;
    timing(8) = peak_value;
    
    threshold = sqrt(threshold*threshold+0.05*0.05); %???what is this ?????

    if threshold > 0.15
        threshold = 0.15;
    end
    if peak_value <= 3*threshold
        if t1_idx > 0
            timing(1) = t1_idx - 1;
            timing(2) = t1_idx - 1;
            timing(3) = t1_idx - 1;
        else
            timing(1) = t1_idx;
            timing(2) = t1_idx;
            timing(3) = t1_idx;
        end
        if t1_idx < last
            timing(5) = t1_idx + 1;
            timing(6) = t1_idx + 1;
            timing(7) = t1_idx + 1;
        else
            timing(5) = t1_idx;
            timing(6) = t1_idx;
            timing(7) = t1_idx;
        end
        return
    end
    
    %t0_idx = FirstLeftIndex(pulse_data_phe,0,t1_idx,1.28*threshold,(3-1.28)*threshold);
    t0_idx = DatPFC_PulseTiming_HeightTiming_Converted_FirstLeftIndex(pulse_data_phe,0,t1_idx,1.28*threshold,(3-1.28)*threshold);
    if t0_idx < 0
        t0_idx = 0;
    end
    timing(1) = t0_idx;
    %t10l_idx = FirstLeftIndex(pulse_data_phe, t0_idx, t1_idx, 0.1*peak_value, threshold);
    t10l_idx = DatPFC_PulseTiming_HeightTiming_Converted_FirstLeftIndex(pulse_data_phe, t0_idx, t1_idx, 0.1*peak_value, threshold);
    timing(2) = t10l_idx;
    %t50l_idx = FirstLeftIndex(pulse_data_phe, t10l_idx, t1_idx, 0.5*peak_value, threshold);
    t50l_idx = DatPFC_PulseTiming_HeightTiming_Converted_FirstLeftIndex(pulse_data_phe, t10l_idx, t1_idx, 0.5*peak_value, threshold);
    timing(3) = t50l_idx;
    
    %t2_idx = FirstRightIndex(pulse_data_phe, last, t1_idx, 1.28*threshold, (3-1.28)*threshold);
    t2_idx = DatPFC_PulseTiming_HeightTiming_Converted_FirstRightIndex(pulse_data_phe, last, t1_idx, 1.28*threshold, (3-1.28)*threshold);
    if t2_idx > last 
        t2_idx = last;
    end
    timing(7) = t2_idx;
    %t10r_idx = FirstRightIndex(pulse_data_phe, t2_idx, t1_idx, 0.1*peak_value, threshold);
    t10r_idx = DatPFC_PulseTiming_HeightTiming_Converted_FirstRightIndex(pulse_data_phe, t2_idx, t1_idx, 0.1*peak_value, threshold);
    timing(6) = t10r_idx;
    %t50r_idx = FirstRightIndex(pulse_data_phe, t10r_idx, t1_idx, 0.5*peak_value, threshold);
    t50r_idx = DatPFC_PulseTiming_HeightTiming_Converted_FirstRightIndex(pulse_data_phe, t10r_idx, t1_idx, 0.5*peak_value, threshold);
    timing(5) = t50r_idx;
    
end