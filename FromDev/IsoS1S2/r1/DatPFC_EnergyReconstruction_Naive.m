function ee = DatPFC_EnergyReconstruction_Naive(ee,XML_Settings,Debug)
if ~exist('Debug','var')
    Debug = 0;
end

%% Initialize variables

    pmt_chs = 1:122;

    % These should all move to iqs in version 2+
    pde = XML_Settings.pde;
    extractioneff = XML_Settings.extractioneff;
    NR_conversion_A = XML_Settings.NR_conversion_A;
    NR_conversion_B = XML_Settings.NR_conversion_B;
    W = XML_Settings.W;
    %{
    [aa b] = size(ee.rqs.event_number);
    ee.rqs.energy_keVee_all = zeros(aa,b);
    ee.rqs.energy_keVee_bot = zeros(aa,b);
    ee.rqs.energy_keVnr_all = zeros(aa,b);
    ee.rqs.energy_keVnr_bot = zeros(aa,b);
    %}

    [c d] = size(ee.rqs.pulse_classification);
    ee.rqs.num_photons        = nan(c,d);
    ee.rqs.num_electrons      = nan(c,d);
    ee.rqs.num_electrons_bot  = nan(c,d);
    ee.rqs.energy_keVee_all   = nan(c,d);
    ee.rqs.energy_keVee_bot   = nan(c,d);
    ee.rqs.energy_keVnr_all   = nan(c,d);
    ee.rqs.energy_keVnr_bot   = nan(c,d);
    ee.rqs.refined_energy_all = nan(c,d);
    ee.rqs.refined_energy_bot = nan(c,d);
    
    if XML_Settings.max_num_pulses == 0
        if Debug
            disp('Fail to run energy reconstruction.')
        end
        ee.info.EnergyReconSuccuss = 0;
        ee.info.EnergyReconError = 'No pulse found.';
        return
    end    
    try
        [dim1, ~] = size(ee.rqs.pulse_classification);

        single_scatter = repmat((sum(ee.rqs.pulse_classification==2).*sum(ee.rqs.pulse_classification==1))==1,[dim1,1]);
        TotalPulseToRecon = ee.rqs.s1s2_pairing(find(ee.rqs.pulse_classification(:,1)==1,1,'last'));
        S1ToUse = zeros(c,d);
        for ii_s1 = 1:TotalPulseToRecon
            S1ToUse( ((ee.rqs.pulse_classification==2) | (ee.rqs.pulse_classification==4)) & (ee.rqs.s1s2_pairing == ii_s1) ) = ...
                ee.rqs.xyz_corrected_pulse_area_all_phe( (ee.rqs.pulse_classification==1) & (ee.rqs.s1s2_pairing == ii_s1) );
        end
        %% IQ fetching (here we go!)
        %{
        single_e_area = [];
        single_e_area_bot = [];
        found_it = 0;

        [a num_iqs] = size(lug_iqs_xml.iq);

         for i=1:num_iqs
          if isfield(lug_iqs_xml.iq(i),'fit')==1   

            if (isfield(lug_iqs_xml.iq(i).fit,'allpmt') && isfield(lug_iqs_xml.iq(i).fit,'botpmt'))==1        
                   single_e_area = lug_iqs_xml.iq(i).fit.allpmt.mean;
                   single_e_area_bot = lug_iqs_xml.iq(i).fit.botpmt.mean;
                   if isfield(lug_iqs_xml.iq(i).fit.botpmt, 'adjrsquare') == 1
                           adjRsqr = lug_iqs_xml.iq(i).fit.botpmt.adjrsquare;
                   else
                           adjRsqr = 1.0;
                   end
                   filename_prefixs = lug_iqs_xml.iq(i).global.filename_prefix;
                   found_it = 1;
                   ee.rqs.energy_found_1eS2_flag = true(aa,b);
                   ee.rqs.energy_bottom_poor_fit = false(aa,b);
            end
          end  
         end
        %} 
         %% Check things are kosher
         %{
         if (found_it == 0)
             % We didn't find our 1eS2 size
             %ee.rqs.energy_found_1eS2_flag = false(aa,b);
             %ee.rqs.energy_bottom_poor_fit = false(aa,b);
             single_e_area = XML_Settings.single_e_area;
             single_e_area_bot = XML_Settings.single_e_area_bot;
        %  elseif ~strcmpi(filename_prefixs,settings.filename_prefix)
        %      % We didn't find our 1eS2 size
        %      ee.rqs.energy_found_1eS2_flag = 0;
        %      ee.rqs.energy_bottom_poor_fit = 0;
        %      single_e_area = XML_Settings.single_e_area;
        %      single_e_area_bot = XML_Settings.single_e_area_bot;
        %  elseif isempty(single_e_area)
        %      % We didn't find our 1eS2 size
        %      ee.rqs.energy_found_1eS2_flag = 0;
        %      ee.rqs.energy_bottom_poor_fit = 0;
        %      single_e_area = XML_Settings.single_e_area;
        %      single_e_area_bot = XML_Settings.single_e_area_bot;
         elseif (adjRsqr < 0.97)
             single_e_area_bot = XML_Settings.single_e_area_bot;
             %ee.rqs.energy_bottom_poor_fit = true(aa,b);
         end
         %}
         single_e_area = XML_Settings.single_e_area;
         single_e_area_bot = XML_Settings.single_e_area_bot;

    %% Calculate Quanta

        ee.rqs.num_photons = (ee.rqs.pulse_classification==1).*ee.rqs.xyz_corrected_pulse_area_all_phe./pde;
        ee.rqs.num_electrons = ((ee.rqs.pulse_classification==2) | (ee.rqs.pulse_classification==4)).*ee.rqs.xyz_corrected_pulse_area_all_phe./(single_e_area.*extractioneff);
        ee.rqs.num_electrons_bot = ((ee.rqs.pulse_classification==2) | (ee.rqs.pulse_classification==4)).*ee.rqs.xyz_corrected_pulse_area_bot_phe./(single_e_area_bot.*extractioneff);

    %% Calculate Energy

        ee.rqs.energy_keVee_all = (ee.rqs.num_electrons + (S1ToUse./pde))*W;
        ee.rqs.energy_keVee_bot = (ee.rqs.num_electrons_bot + (S1ToUse./pde))*W;
        ee.rqs.energy_keVnr_all = (ee.rqs.energy_keVee_all/NR_conversion_A).^(1/NR_conversion_B);
        ee.rqs.energy_keVnr_bot = (ee.rqs.energy_keVee_bot/NR_conversion_A).^(1/NR_conversion_B);
        
    %% Refine Energy (keVee)

        %ee.rqs.refined_energy_all = (repmat(Beta(ee.rqs.energy_keVee_all),dim1,1).*ee.rqs.num_photons + ee.rqs.num_electrons).*W.*Alpha(ee.rqs.energy_keVee_all);
        %ee.rqs.refined_energy_bot = (repmat(Beta(ee.rqs.energy_keVee_bot),dim1,1).*ee.rqs.num_photons + ee.rqs.num_electrons_bot).*W.*Alpha(ee.rqs.energy_keVee_bot);
        ee.rqs.refined_energy_all = (Beta(ee.rqs.energy_keVee_all).*(S1ToUse./pde) + ee.rqs.num_electrons).*W.*Alpha(ee.rqs.energy_keVee_all);
        ee.rqs.refined_energy_bot = (Beta(ee.rqs.energy_keVee_bot).*(S1ToUse./pde) + ee.rqs.num_electrons_bot).*W.*Alpha(ee.rqs.energy_keVee_bot);
        
        ee.info.EnergyReconSuccuss = 1;
        ee.info.EnergyReconError = '';
    catch exception
        if Debug
            disp('Fail to run energy reconstruction.')
        end
        ee.info.EnergyReconSuccuss = 0;
        ee.info.EnergyReconError = exception.identifier;
    end
end

function betas = Beta(raw_energy_all)
    betas = 0.83827*ones(size(raw_energy_all)) - 0.83336 * exp(-0.46154 * raw_energy_all);
end

function alphas = Alpha(raw_energy_all)
    alphas =1.0504*ones(size(raw_energy_all)) + 1.0054*Beta(raw_energy_all) ...
        - 1.8906*(Beta(raw_energy_all)).^2 + 1.277*(Beta(raw_energy_all)).^3 ...
        - 0.41817*(Beta(raw_energy_all)).^4;
end



