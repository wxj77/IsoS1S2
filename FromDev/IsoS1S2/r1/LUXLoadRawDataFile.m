function [raw livetime xlm settings] = LUXLoadRawDataFile( filename, data_path , options)
%
%function [raw livetime xlm settings] = LUXLoadRawDataFile( filename, data_path )
%
%This function loads the raw .dat file into memory. This is useful if your acquisition is in multi-event
%mode and you want to look at everything the Struck saves, ignoring event triggers
%
% inputs
%
%  filename 
%  data_path 
%  options.no_warnings - suppress print out of warnings. Errors still report.
%
% e.g. clear options; options.no_warnings = 1;
%      [raw livetime xlm] = LUXLoadRawDataFile( sprintf('%s_f%09d.dat',series,ii) , dat_path , options);
%
%
% 2012-02-13 JJC - now gunzips and rezips .gz files
% 2012-11-19 JJC - changed output from pulse to pod.
% 2012-12-25 MM - change format when read out XLM hit vectors from int8 to uint8 
% 2013-01-13 RJG - add options. Allow suppression of gunzip warning
% notifcation. Helps when iterating over many files, and don't want to see
% warnings
% 2013-07-11 RJG - Minor bug fix on warnings switch
% 2014-06-13 MM - Now read out xlm_trigseqnum and ccheck correctly in line 124 and 133.
%

 if nargin<3
     options.no_warnings=0;
 end
 if nargin==3
    if ~isfield(options, 'no_warnings')
        options.no_warnings=0;
    end
 end
 
raw = [];
livetime = [];
xlm = [];

if ~exist('filename', 'var') || isempty(filename)
	fprintf('you must enter a filename\n');
	return;
end


if ~exist('data_path', 'var')
	data_path = './';
    fprintf('you did not enter a data_path, assuming ./\n');
end

file_list = dir(data_path);
file_found = 0;
for ii = 1:length(file_list)
    if strcmp(filename, file_list(ii).name)
    file_found = 1;
    gz_file = 0;
    end
end
if file_found==0
    if ~options.no_warnings, fprintf('did not find file %s in %s\nchecking for gzipped file...\n',filename,data_path), end
    for ii = 1:length(file_list)
        if strcmp([filename '.gz'], file_list(ii).name)
            file_found = 1;
            gz_file = 1;
        end
    end
end

if file_found==0
    fprintf('did not find file %s (gz or otherwise) in %s\n',filename,data_path);
    return;
end


	

	% Begin Reading File %
	
        filename_with_path = [data_path,'/',filename];
        if gz_file==1
           copyfile([filename_with_path '.gz'], './'); 
           fprintf('gunzipping %s ...\n', [filename_with_path '.gz']);
           gunzip(['./' filename '.gz']);
           fid = fopen(['./' filename],'rb','l');
        else
           fid = fopen(filename_with_path,'rb','l');
        end
        
        if fid<0
            fprintf('Failed to find the file (%s)\n',filename_with_path);
            return
        end

        
        % Read dat xml file
        endianness = fread(fid,1,'ulong');
        xml_length = fread(fid,1,'ulong');
        xml_string = char(fread(fid,xml_length,'char'))';
        settings = XMLParser(char(xml_string));
        
        
        % / Read dat xml file
        
		%endianness = fread(fid, 1, 'uint32');
		% read file header
		file_header = fread(fid,8,'uint32');
		nb_chs = file_header(6);
		%if ~exist('nb_seqs','var') % if first file, initialize variables
		%	for ch=1:nb_chs
		%		pulse_data = [];
		%		baseline = [];
		%		timestamp = [];
		%	end
		%end
		if isfield(settings, 'daq_settings')
        if isfield(settings.daq_settings.sis3301.global, 'read_xlm')
        if settings.daq_settings.sis3301.global.read_xlm==1
            xlm.xlm_nb_triggers = fread(fid,1,'uint16');
            xlm.xlm_nb_ddcs = fread(fid,1,'uint16');
            for n_xlm=1:xlm.xlm_nb_triggers
               xlm.xlm_trig_timestamp(n_xlm) = fread(fid,1,'uint64');
               if settings.daq_settings.global.daq_version>7.0
                  xlm.xlm_trigseqnum(n_xlm) = fread(fid,1,'uint32'); 
               end
               xlm.xlm_max_filter_response(n_xlm) = fread(fid,1,'uint32');
               xlm.xlm_max_ch_ID(n_xlm) = fread(fid,1,'char');
               for d_xlm=1:xlm.xlm_nb_ddcs
                  xlm.xlm_S1_hit_vector(n_xlm,d_xlm) = fread(fid,1,'uint8');
                  xlm.xlm_S2_hit_vector(n_xlm,d_xlm) = fread(fid,1,'uint8');
               end
               if settings.daq_settings.global.daq_version>7.0
                  xlm.xlm_ccheck(n_xlm) = fread(fid,1,'uint8'); 
               end
            end
        end
        end
        end
		% read livetime header
		nb_seqs = fread(fid,1,'uint16');
		%keyboard
        %livetime_latch = zeros(1,nb_seqs);
		%livetime_end = zeros(1,nb_seqs);
		for iseq=1:nb_seqs
			livetime(iseq).latch = fread(fid,1,'uint64');
			livetime(iseq).end = fread(fid,1,'uint64');
		end
		
		% loop over channels
		for ch=1:nb_chs
                
				ps_tot(ch) = 0;
			
			channel_header = fread(fid,4,'uint16');
			nb_pulses = channel_header(2);
			
			% loop over pulses
			for ps=1:nb_pulses
			
				pulse_header = fread(fid,9,'uint32');
				pulse_length = pulse_header(7);
                pulse_trigger_flags = pulse_header(8);
				pulse_baseline = pulse_header(9);
				pulse_timestamp = pulse_header(6) + pulse_header(5)*2^32;
				pulse_samples = fread(fid,pulse_length,'uint16')';
				%pulse_samples=pulse_samples(2:end); %kludge to fix first sample from acq
				%pulse_length=pulse_length-1;
				
				%{
				pulse_baseline_vec = ones(1,pulse_length)*pulse_baseline;
				pulse_timestamp_vec = ones(1,pulse_length)*(pulse_timestamp-1) + (1:1:pulse_length);
				
				if file_number==file_numbers_list(1) && ps==1 % first file
					raw.ch(ch).pulse_data = pulse_samples; % add new pulse on to end of previous pulses
					raw.ch(ch).baseline = pulse_baseline_vec;
					raw.ch(ch).timestamp = pulse_timestamp_vec;
				else
					raw.ch(ch).pulse_data = [raw.ch(ch).pulse_data pulse_samples]; % add new pulse on to end of previous pulses
					raw.ch(ch).baseline = [raw.ch(ch).baseline pulse_baseline_vec];
					raw.ch(ch).timestamp = [raw.ch(ch).timestamp pulse_timestamp_vec];
				end
				%}
				ps_tot(ch) = ps_tot(ch) + 1;
				raw.ch(ch).pod(ps_tot(ch)).pod_data = pulse_samples';
				%raw.ch(ch).pulse(ps_tot(ch)).pulse_data_mV = pulse_samples'.*(2000/2^14);
				raw.ch(ch).pod(ps_tot(ch)).baseline = pulse_baseline;
				raw.ch(ch).pod(ps_tot(ch)).baseline_mV = pulse_baseline.*(2000/2^14);
				raw.ch(ch).pod(ps_tot(ch)).pod_data_mV = pulse_samples'.*(2000/2^14) - pulse_baseline.*(2000/2^14);
                raw.ch(ch).pod(ps_tot(ch)).timestamp = pulse_timestamp;
				raw.ch(ch).pod(ps_tot(ch)).length = pulse_length;
                raw.ch(ch).pod(ps_tot(ch)).trigger_flags = pulse_trigger_flags;
                raw.ch(ch).pod(ps_tot(ch)).start = pulse_timestamp - livetime(iseq).latch;
			end
			
		end
		
		fclose(fid);
        if gz_file==1
            delete(['./' filename]);
            delete(['./' filename '.gz']);
            %fprintf('done with %s, gzipping it back...\n', filename_with_path);
            %gzip([filename_with_path]);
        end

