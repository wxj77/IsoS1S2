function ee = DatPFC_PosRecon(ee,XML_Settings,Iq,SR_Version,Debug)
if ~exist('SR_Version','var')
    disp('SR_Version is not specified. Use default SR2.0.')
    SR_Version = 'SR2.0';
end
if ~exist('Debug','var')
    Debug = 0;
end

try
    if strcmp(SR_Version,'SR2.0')
        
        rec_set = XML_Settings;

        %% Input checking for the structure rec_set.
        %%

        if ~isfield(rec_set, 'file_lrfmat'),
            %positions = findstr(data_path_evt, '/');    
            %rec_set.file_lrfmat = [char(data_path_evt(1:positions(end))) '.' filename_evt(1:end-20) '_rec_fun.mat'];%point to a specific _rec_fun.mat
            %rec_set.file_lrfmat = ['G:\LUXData\dat\lux10_20130728T2255\rec_fun.mat'];
        end

        rec_set.PMTS_Working = rec_set.PMTS_To_Use==1; 
        rec_set.search_for_double_scatterer = 0; 
        rec_set.compute_map = 0; 
        rec_set.tests = 0; 
        rec_set.only_radial_component = 0; 
        rec_set.get_estimated_peak_areas = 1;

        %rec_set.saturated_pmt_phe_limit
        rec_set.saturated_pmt_phe_limit = [10000  12000  10000   7000   8000   7000  12000  10000   2500   9000   7000   8000  12000   7000   6000   6000   6000   6000   6000   8000  10000   8000 12000   6000   8000  12000   7000   8000   7000   6000  12000   8000   6000  12000  12000   7000   8000   7000   5000   5500  12000  12000   7000   8000   7000  11000  12000  10000  12000   6000   6000   4000   5000  12000   7000  12000  12000  12000  12000   4000  12000];


        %% Prepare LRFs

        id_prog = 'PositionReconstruction_MercuryI: '; 

        %{
        for ii = 1:length(Iq.lrfs.iq)
            if isfield(Iq.lrfs.iq(ii).global, 'iq_type') 
                if strcmp(Iq.lrfs.iq(ii).global.iq_type,'lrfs') & isfield(Iq.lrfs.iq(ii).global, 'algorithm_name')
                    if strcmp(Iq.lrfs.iq(ii).global.algorithm_name,'functional') == 1 | strcmp(Iq.lrfs.iq(ii).global.algorithm_name,rec_set.algorithm) 
                        lrf_iq = Iq.lrfs.iq(ii);
                    end
                end
            end
        end
        %}
        if strcmp(Iq.lrfs.iq.global.algorithm_name,'functional') == 1 || strcmp(Iq.lrfs.iq.global.algorithm_name,rec_set.algorithm) 
            lrf_iq = Iq.lrfs.iq;
        end


        %% Get PMT Gains
        if rec_set.MLM_maxphe > 0
            for ii = 1:length(Iq.pmt_gains.iq)
                if isfield(Iq.pmt_gains.iq(ii).global, 'iq_type') 
                    if strcmp(Iq.pmt_gains.iq(ii).global.iq_type,'pmt_gains')
                        lrf_iq.sphe = [Iq.pmt_gains.iq(ii).fit.channel(1:122).sigma_mVns_per_phe]./[Iq.pmt_gains.iq(ii).fit.channel(1:122).mVns_per_phe];
                        lrf_iq.sphe(~inrange(lrf_iq.sphe, 0.1, 2.)) = 0.6;
                    end
                end
            end
            if ~isfield(lrf_iq, 'sphe')
                %fprintf(['PositionReconstruction_MercuryI. We did not found the PMT gains in the IQ settings file.\n']);
                %fprintf(['The chi2 method will be used instead.\n']);  
                rec_set.MLM_maxphe = -100000;
            end
        end
        %% Selection of the events for the reconstruction. Previously in an
        %% independent file

        NPULS = size(ee.rqs.peak_area_phe, 1);
        if NPULS == 0
            ee.rqs.reconstructed  = nan(NPULS, 1);
            ee.rqs.rec_energy_flag = zeros(NPULS, 1);

            ee.rqs.x_cm = 100*ones(size(ee.rqs.pulse_area_phe)); % a very out of range value
            ee.rqs.y_cm = 100*ones(size(ee.rqs.pulse_area_phe)); % a very out of range value
            ee.rqs.minimization_flag = -1*ones(size(ee.rqs.pulse_area_phe)); % a failure value
            ee.rqs.chi2 = -1*ones(size(ee.rqs.pulse_area_phe)); % a failure value


            if Debug
                disp('Fail to run position reconstruction.')
            end

            ee.info.PosReconSuccess = 0;
            ee.info.PosReconError = 'No pulse found';
            ee.info.PosReconVersion = SR_Version;
            return
        end
        ee.rqs.reconstructed  = ones(NPULS, 1);

        NumberOfTopPMTsWithSignal = squeeze(sum(ee.rqs.peak_area_phe(:, [1:60 121],:)>0, 2));
        ee.rqs.reconstructed(squeeze(sum(ee.rqs.peak_area_phe,2))<rec_set.min_num_of_phe) = 0;
        ee.rqs.rec_energy_flag = zeros(NPULS, 1);

        if rec_set.energy_minimization
            ee.rqs.rec_energy_flag = squeeze(sum(ee.rqs.peak_area_phe,2))>rec_set.energy_minimum_for_minimization;
        else
            ee.rqs.rec_energy_flag = ee.rqs.reconstructed.*0;
        end

        %%% 
        %%% Array equalization and the preliminary estimative of the position of the event the position of the event.
        %%%

        if ~isfield(lrf_iq, 'sigma_zero') 
            lrf_iq.sigma_zero = 0.4; 
        end

        if isfield(ee.rqs, 'pulse_area_phe')
            [ee.rqs rec_set] = DatPFC_MercuryPrepareMinimization(ee.rqs, rec_set, lrf_iq);
            %ee.rqs.reconstructed = double(ee.rqs.reconstructed);
            ee.rqs = DatPFC_PositionReconstruction_MercuryGetXY(ee.rqs, rec_set, lrf_iq);
            ee.rqs = rmfield(ee.rqs, {'x_cm_old', 'y_cm_old', 'minimization_flag', 'rec_energy_flag', 'peak_area_eq', 'pulse_area_eq'});
        else
            %disp(fprintf('\n **** %s 0 Events in the rq file %s. **** \n',myname, filename_rq));
            ee.rqs.x_cm = 100*ones(size(ee.rqs.pulse_area_phe)); % a very out of range value
            ee.rqs.y_cm = 100*ones(size(ee.rqs.pulse_area_phe)); % a very out of range value
            ee.rqs.minimization_flag = -1*ones(size(ee.rqs.pulse_area_phe)); % a failure value
            ee.rqs.chi2 = -1*ones(size(ee.rqs.pulse_area_phe)); % a failure value
        end

        if rec_set.get_estimated_peak_areas==0
            ee.rqs = rmfield(ee.rqs, 'peak_area_rec');
        end

        ee.info.PosReconSuccess = 1;
    	ee.info.PosReconError = '';
    	ee.info.PosReconVersion = SR_Version;
    elseif strcmp(SR_Version,'SR2.1')
       
        rec_set = XML_Settings;

        %% Input checking for the structure rec_set.
        %%

        if ~isfield(rec_set, 'file_lrfmat'),
            %positions = findstr(data_path_evt, '/');    
            %rec_set.file_lrfmat = [char(data_path_evt(1:positions(end))) '.' filename_evt(1:end-20) '_rec_fun.mat'];%point to a specific _rec_fun.mat
            %rec_set.file_lrfmat = ['G:\LUXData\dat\lux10_20130728T2255\rec_fun.mat'];
        end


        %%
        if ~isfield(rec_set, 'PMTsOff')
            rec_set.PMTsOff =  XML_Settings.pmts_off; 
        end
        
        rec_set.PMTS_To_Use = ones(1, 122);
        rec_set.PMTS_To_Use(rec_set.PMTsOff) = 0;


        rec_set.PMTS_Working = rec_set.PMTS_To_Use==1; rec_set.search_for_double_scatterer = 0; rec_set.compute_map = 0; rec_set.tests = 0; rec_set.only_radial_component = 0; rec_set.get_estimated_peak_areas = 1;

        %rec_set.saturated_pmt_phe_limit
        rec_set.saturated_pmt_phe_limit = [10000  12000  10000   7000   8000   7000  12000  10000   2500   9000   7000   8000  12000   7000   6000   6000   6000   6000   6000   8000  10000   8000 12000   6000   8000  12000   7000   8000   7000   6000  12000   8000   6000  12000  12000   7000   8000   7000   5000   5500  12000  12000   7000   8000   7000  11000  12000  10000  12000   6000   6000   4000   5000  12000   7000  12000  12000  12000  12000   4000  12000];


        %% Prepare LRFs

        id_prog = 'PositionReconstruction_MercuryI: '; 

        %{
        for ii = 1:length(Iq.lrfs.iq)
            if isfield(Iq.lrfs.iq(ii).global, 'iq_type') 
                if strcmp(Iq.lrfs.iq(ii).global.iq_type,'lrfs') & isfield(Iq.lrfs.iq(ii).global, 'algorithm_name')
                    if strcmp(Iq.lrfs.iq(ii).global.algorithm_name,'functional') == 1 | strcmp(Iq.lrfs.iq(ii).global.algorithm_name,rec_set.algorithm) 
                        lrf_iq = Iq.lrfs.iq(ii);
                    end
                end
            end
        end
        %}
        if strcmp(Iq.lrfs.iq.global.algorithm_name,'functional') == 1 || strcmp(Iq.lrfs.iq.global.algorithm_name,rec_set.algorithm) 
            lrf_iq = Iq.lrfs.iq;
        end


        %% Get PMT Gains
        if rec_set.MLM_maxphe > 0
            for ii = 1:length(Iq.pmt_gains.iq)
                if isfield(Iq.pmt_gains.iq(ii).global, 'iq_type') 
                    if strcmp(Iq.pmt_gains.iq(ii).global.iq_type,'pmt_gains')
                        lrf_iq.sphe = [Iq.pmt_gains.iq(ii).fit.channel(1:122).sigma_mVns_per_phe]./[Iq.pmt_gains.iq(ii).fit.channel(1:122).mVns_per_phe];
                        lrf_iq.sphe(~inrange(lrf_iq.sphe, 0.1, 2.)) = 0.6;
                    end
                end
            end
            if ~isfield(lrf_iq, 'sphe')
                %fprintf(['PositionReconstruction_MercuryI. We did not found the PMT gains in the IQ settings file.\n']);
                %fprintf(['The chi2 method will be used instead.\n']);  
                rec_set.MLM_maxphe = -100000;
            end
        end
        %% Selection of the events for the reconstruction. Previously in an
        %% independent file

        NPULS = size(ee.rqs.peak_area_phe, 1);
        if NPULS == 0
            ee.rqs.reconstructed  = nan(NPULS, 1);
            ee.rqs.rec_energy_flag = zeros(NPULS, 1);

            ee.rqs.x_cm = 100*ones(size(ee.rqs.pulse_area_phe)); % a very out of range value
            ee.rqs.y_cm = 100*ones(size(ee.rqs.pulse_area_phe)); % a very out of range value
            ee.rqs.minimization_flag = -1*ones(size(ee.rqs.pulse_area_phe)); % a failure value
            ee.rqs.chi2 = -1*ones(size(ee.rqs.pulse_area_phe)); % a failure value


            if Debug
                disp('Fail to run position reconstruction.')
            end

            ee.info.PosReconSuccess = 0;
            ee.info.PosReconError = 'No pulse found';
            ee.info.PosReconVersion = SR_Version;
            return
        end
        
        ee.rqs.reconstructed  = ones(NPULS, 1);

        NumberOfTopPMTsWithSignal = squeeze(sum(ee.rqs.peak_area_phe(:, [1:60 121],:)>0, 2));
        ee.rqs.reconstructed(squeeze(sum(ee.rqs.peak_area_phe,2))<rec_set.min_num_of_phe) = 0;
        ee.rqs.rec_energy_flag = zeros(NPULS, 1);

        if rec_set.energy_minimization
            ee.rqs.rec_energy_flag = squeeze(sum(ee.rqs.peak_area_phe,2))>rec_set.energy_minimum_for_minimization;
        else
            ee.rqs.rec_energy_flag = ee.rqs.reconstructed.*0;
        end

        %%% 
        %%% Array equalization and the preliminary estimative of the position of the event the position of the event.
        %%%

        if ~isfield(lrf_iq, 'sigma_zero') 
            lrf_iq.sigma_zero = 0.4; 
        end

        if isfield(ee.rqs, 'pulse_area_phe')
            [ee.rqs rec_set] = DatPFC_MercuryPrepareMinimization(ee.rqs, rec_set, lrf_iq);
            %ee.rqs.reconstructed = double(ee.rqs.reconstructed);
            ee.rqs = DatPFC_PositionReconstruction_MercuryGetXY(ee.rqs, rec_set, lrf_iq);
            ee.rqs = rmfield(ee.rqs, {'x_cm_old', 'y_cm_old', 'minimization_flag', 'rec_energy_flag', 'peak_area_eq', 'pulse_area_eq'});
        else
            %disp(fprintf('\n **** %s 0 Events in the rq file %s. **** \n',myname, filename_rq));
            ee.rqs.x_cm = 100*ones(size(ee.rqs.pulse_area_phe)); % a very out of range value
            ee.rqs.y_cm = 100*ones(size(ee.rqs.pulse_area_phe)); % a very out of range value
            ee.rqs.minimization_flag = -1*ones(size(ee.rqs.pulse_area_phe)); % a failure value
            ee.rqs.chi2 = -1*ones(size(ee.rqs.pulse_area_phe)); % a failure value
        end

        if rec_set.get_estimated_peak_areas==0
            ee.rqs = rmfield(ee.rqs, 'peak_area_rec');
        end

        ee.info.PosReconSuccess = 1;
        ee.info.PosReconError = '';
    	ee.info.PosReconVersion = SR_Version;
    else
        error('Do not recognize SR version.');
    end
catch exception
    NPULS = size(ee.rqs.peak_area_phe, 1);

    ee.rqs.reconstructed  = nan(NPULS, 1);
    ee.rqs.rec_energy_flag = zeros(NPULS, 1);

    ee.rqs.x_cm = 100*ones(size(ee.rqs.pulse_area_phe)); % a very out of range value
    ee.rqs.y_cm = 100*ones(size(ee.rqs.pulse_area_phe)); % a very out of range value
    ee.rqs.minimization_flag = -1*ones(size(ee.rqs.pulse_area_phe)); % a failure value
    ee.rqs.chi2 = -1*ones(size(ee.rqs.pulse_area_phe)); % a failure value
    
    
    if Debug
        disp('Fail to run position reconstruction.')
    end

    ee.info.PosReconSuccess = 0;
    ee.info.PosReconError = exception.identifier;
    ee.info.PosReconVersion = SR_Version;
end
