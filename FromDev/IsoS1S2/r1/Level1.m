function [] = Level1(DataSets,RunSetting)

	disp(size(DataSets,2))
	parfor dn = 1:size(DataSets,2)
		DataSet=DataSets{1,dn};
		cd '/home/dkhaitan/IsolatedS1S2/';
		DataPath = ['/home/dkhaitan/LUXData/IsolatedS1S2/' DataSet '/'];
		DataOutPath = ['/home/dkhaitan/IsolatedS1S2/test_out/' DataSet '/'];
		MaxDriftLength = 35000; 
		FileList = dir([DataPath '*.dat']);
		FileListToRun = length(FileList);


		if ~exist(strrep(DataPath,'dat','DatPFC_Extended'),'dir')
		    mkdir(strrep(DataPath,'dat','DatPFC_Extended'));
		end

		Setting = XMLReader_framework(RunSetting);
		isolated_s1_rate = 0;
		isolated_s2_rate = 0;
		s1_spect = zeros(1,1000);
		s2_spect = zeros(1,1000);
		livetime_total=0;


		%Setting = xmlread(RunSetting);
		for ii = 1:1000
			try
				disp(FileList(ii).name);
				[FileOut,ee_merge] = Level2(FileList(ii).name,DataPath,DataOutPath,DataSet,Setting);
				[temp_isolated_s1_rate,temp_isolated_s2_rate,temp_s1_spect,temp_s2_spect,livetime] = IsolateS1S2s(ee_merge,MaxDriftLength);

				isolated_s1_rate = (livetime_total*isolated_s1_rate + livetime*temp_isolated_s1_rate)/(livetime_total + livetime);
				isolated_s2_rate = (livetime_total*isolated_s2_rate + livetime*temp_isolated_s2_rate)/(livetime_total + livetime);
				s1_spect = (s1_spect*livetime_total + livetime*temp_s1_spect)/(livetime_total + livetime);
				s2_spect = (s2_spect*livetime_total + livetime*temp_s2_spect)/(livetime_total + livetime);
				livetime_total=livetime_total + livetime;
				disp([temp_isolated_s1_rate,temp_isolated_s2_rate,livetime]);

			catch exception
				disp(['Fail ' FileList(ii).name])
			end
			clearvars -global -except DataSets DataSet DataPath DataOutPath MaxDriftLength FileList Setting isolated_s1_rate isolated_s2_rate s1_spect s2_spect livetime 
		end
		fid_s2_spect=fopen([DataOutPath 's2_spect.csv'],'at');
		fprintf(fid_s2_spect,'%s,',DataSet);
		fprintf(fid_s2_spect,'%s,',num2str(isolated_s2_rate));
		dlmwrite([DataOutPath 's2_spect.csv'],s2_spect,'-append','delimiter',',');
		fclose(fid_s2_spect);

		fid_s1_spect=fopen([DataOutPath 's1_spect.csv'],'at');
		fprintf(fid_s1_spect,'%s,',DataSet);
		fprintf(fid_s1_spect,'%s,',num2str(isolated_s1_rate));
		dlmwrite([DataOutPath 's1_spect.csv'],s1_spect,'-append','delimiter',',');
		fclose(fid_s1_spect);

	end
end
