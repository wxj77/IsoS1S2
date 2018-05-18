function left = DatPFC_PulseTiming_HeightTiming_Converted_FirstLeftIndex(pulse_data_phe,start,stop,above,noise)
% "start" and "stop" are cpp index.
%
    if start < (stop-1)
        left0 =  find(...
            (...
            (pulse_data_phe( (start:(stop-2))+1+1 ) > above )...% the +1 is from looking at the next entry the +1 is for shifting index
            &...
            (pulse_data_phe( (start:(stop-2))+2+1 ) > above )...% the +2 is from looking at the second next entry the +1 is for shifting index
            )...
            |...
            (pulse_data_phe( (start:(stop-2))+1+1 ) >( above+3*noise) )...% the +1 is from looking at the next entry the +1 is for shifting index
            ,1,'first');  
        if ~isempty(left0)
            Temp = start:(stop-2);
            left = Temp(left0);%We don't need the -1 is to shift the index back to c++ convention here, right ?  
        else
            left = stop-1;
        end
    else
        left = start;
    end
end