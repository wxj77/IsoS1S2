function right = DatPFC_PulseTiming_HeightTiming_Converted_FirstRightIndex(pulse_data_phe,start,stop,above,noise)
    if start > (stop+1)
        right0 = find(... 
            (...
            (pulse_data_phe( ((stop+2):start)-1+1 ) > above )...% the -1 is from looking at the next entry the +1 is for shifting index
            &...
            (pulse_data_phe( ((stop+2):start)-2+1 ) > above )...% the -2 is from looking at the second next entry the +1 is for shifting index
            )...
            |...
            (pulse_data_phe( ((stop+2):start)-1+1 ) > (above+3*noise) )...% the -1 is from looking at the next entry the +1 is for shifting index
            ,1,'last');

        if ~isempty(right0)
            Temp = (stop+2):start;
            right = Temp(right0);%We don't need the -1 is to shift the index back to c++ convention here, right ?  
        else
            right = stop+1;
            %right = start;
        end
    else
        right = start;
        %right = stop+2;
    end
end
