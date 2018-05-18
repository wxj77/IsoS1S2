function ee = DatPFC_Corrections_PositionCorrection(ee,XML_Settings,Iq,SR_Version,Debug)
if ~exist('SR_Version','var')
    disp('SR_Version is not specified. Use default SR2.0.')
    SR_Version = 'SR2.0';
end 
if ~exist('Debug','var')
    Debug = 0;
end

try
    lrf_iq = Iq.xy_rec_cor.iq;
    
    if ee.rqs.num_pulses_found == 0
        %This one is the original complaint.
        %disp(fprintf('\n\n %s: The IQ xy_rec_cor was not found\n\n************* Fatal Error*************\n\n',myname));
        ee.rqs.x_corrected = nan(size(ee.rqs.pulse_area_phe));
        ee.rqs.y_corrected = nan(size(ee.rqs.pulse_area_phe));
        if Debug
            disp('Fail to run position correction.')
        end
        ee.info.PosCorrSuccess = 0;
        ee.info.PosCorrError = 'No pulse found.';
        ee.info.PosCorrVersion = SR_Version;
        return
    end

    if ~isstruct(lrf_iq)
        %This one is the original complaint.
        %disp(fprintf('\n\n %s: The IQ xy_rec_cor was not found\n\n************* Fatal Error*************\n\n',myname));
        ee.rqs.x_corrected = nan(size(ee.rqs.pulse_area_phe));
        ee.rqs.y_corrected = nan(size(ee.rqs.pulse_area_phe));
        if Debug
            disp('Fail to run position correction.')
        end
        ee.info.PosCorrSuccess = 0;
        ee.info.PosCorrError = 'The IQ xy_rec_cor was not found.';
        ee.info.PosCorrVersion = SR_Version;
        return
    end

    if ~isfield(lrf_iq, 'file_table_Mercury') && isfield(lrf_iq, 'file_table')
        lrf_iq.file_table_Mercury = lrf_iq.file_table;
    end

    if isfield(lrf_iq, 'file_table_Mercury')
        if ~any(lrf_iq.file_table_Mercury)
            %This one is the original complaint.
            %disp(['Can not find ' lrf_iq.file_table_Mercury]);
            ee.rqs.x_corrected = nan(size(ee.rqs.pulse_area_phe));
            ee.rqs.y_corrected = nan(size(ee.rqs.pulse_area_phe));
            if Debug
                disp('Fail to run position correction.')
            end
            ee.info.PosCorrSuccess = 0;
            ee.info.PosCorrError = ['Can not find ' lrf_iq.file_table_Mercury];
            ee.info.PosCorrVersion = SR_Version;
        else
            table_corrections = load(lrf_iq.file_table_Mercury);
            %% Run the function. 
            ee.rqs = DatPFC_Corrections_PositionCorrection_Function(ee.rqs, table_corrections,SR_Version);
            ee.info.PosCorrSuccess = 1;
            ee.info.PosCorrError = '';
            ee.info.PosCorrVersion = SR_Version;
        end
    else
        ee.rqs.x_corrected = nan(size(ee.rqs.pulse_area_phe));
        ee.rqs.y_corrected = nan(size(ee.rqs.pulse_area_phe));
        if Debug
            disp('Fail to run position correction.')
        end
        ee.info.PosCorrSuccess = 0;
        ee.info.PosCorrError = 'Can not find calibration table.';
        ee.info.PosCorrVersion = SR_Version;
    end
    
catch exception
    ee.rqs.x_corrected = nan(size(ee.rqs.pulse_area_phe));
    ee.rqs.y_corrected = nan(size(ee.rqs.pulse_area_phe));
    if Debug
        disp('Fail to run position correction.')
    end
    ee.info.PosCorrSuccess = 0;
    ee.info.PosCorrError = exception.identifier;
    ee.info.PosCorrVersion = SR_Version;
end
%{
if isfield(lrf_iq, 'file_table_TAXY')
    position_correction_table_path = [position_correction_path(1:(end-numel(myname)-2)) lrf_iq.file_table_TAXY];
    if ~any(position_correction_table_path)
        disp(fprintf('\n\n %s: The file %s was not found in %s\n\n*************\n\n',myname, lrf_iq.file_table, position_correction_path(1:(end-numel(myname)-2))));
    else
        table_corrections = load(position_correction_table_path);
        %% Run the function. 
        dp = Corrections_PositionCorrection_Function(dp, table_corrections);        
    end
    dp.x_tmplt_corrected = single(dp.x_tmplt_corrected);
    dp.y_tmplt_corrected = single(dp.y_tmplt_corrected);
    dp.xy_sigma_corrected = single(dp.xy_sigma_corrected);
end
%}

