
function Level1_sbatch_2(iij)

iij = 1402+iij
pc=parcluster('local');
JOB_ID=getenv('SLURM_JOBID');
CPUS=str2num(getenv('SLURM_JOB_CPUS_PER_NODE'));
pc.JobStorageLocation=strcat('/local_scratch/',JOB_ID);
parpool(pc,CPUS)

%{
DataSets{1,1}='lux10_20140903T1918';
DataSets{1,2}='lux10_20140903T2300';
DataSets{1,3}='lux10_20140903T1918_2';
DataSets{1,4}='lux10_20140903T2300_2';
%%
DataSets{1,1}='lux10_20140904T0325';
DataSets{1,2}='lux10_20140911T1548';
DataSets{1,3}='lux10_20140911T2113';
DataSets{1,4}='lux10_20140912T0247';
DataSets{1,5}='lux10_20140912T0822';
DataSets{1,6}='lux10_20140912T1357';
DataSets{1,7}='lux10_20140912T1922';
DataSets{1,8}='lux10_20140913T0047';
DataSets{1,9}='lux10_20140913T0607';
DataSets{1,10}='lux10_20141127T1910';
DataSets{1,11}='lux10_20141127T2303';
DataSets{1,12}='lux10_20141128T0257';
DataSets{1,13}='lux10_20141128T0656';
DataSets{1,14}='lux10_20141128T1030';
DataSets{1,15}='lux10_20141129T1514';
DataSets{1,16}='lux10_20141129T1838';
DataSets{1,17}='lux10_20141129T2207';
DataSets{1,18}='lux10_20141130T0051';
DataSets{1,19}='lux10_20141130T0531';
DataSets{1,20}='lux10_20141130T0908';
DataSets{1,21}='lux10_20141130T1250';
DataSets{1,22}='lux10_20141130T1635';
DataSets{1,23}='lux10_20141130T2023';
DataSets{1,24}='lux10_20141201T0009';
DataSets{1,25}='lux10_20141201T0354';
DataSets{1,26}='lux10_20141201T0927';
DataSets{1,27}='lux10_20141201T1305';
DataSets{1,28}='lux10_20141201T1701';
DataSets{1,29}='lux10_20141201T2054';
DataSets{1,30}='lux10_20141202T0048';
DataSets{1,31}='lux10_20141202T0442';
DataSets{1,32}='lux10_20141202T0833';
DataSets{1,33}='lux10_20141202T1129';
DataSets{1,34}='lux10_20141202T1449';
DataSets{1,35}='lux10_20141202T1835';
DataSets{1,36}='lux10_20141202T2245';
DataSets{1,37}='lux10_20141203T0310';
DataSets{1,38}='lux10_20141203T2205';
DataSets{1,39}='lux10_20141204T0151';
DataSets{1,40}='lux10_20141204T0541';
DataSets{1,41}='lux10_20141204T0913';
DataSets{1,42}='lux10_20141204T1329';
DataSets{1,43}='lux10_20141204T1720';
DataSets{1,44}='lux10_20141204T2205';
DataSets{1,45}='lux10_20141205T0254';
DataSets{1,46}='lux10_20141205T0751';
DataSets{1,47}='lux10_20141205T1239';
DataSets{1,48}='lux10_20141205T1636';
DataSets{1,49}='lux10_20141205T2135';
DataSets{1,50}='lux10_20141206T0232';
DataSets{1,51}='lux10_20141206T0732';
DataSets{1,52}='lux10_20141206T1233';
DataSets{1,53}='lux10_20141206T1710';
DataSets{1,54}='lux10_20141206T2206';
DataSets{1,55}='lux10_20141207T0051';
DataSets{1,56}='lux10_20141207T0548';
DataSets{1,57}='lux10_20141207T1046';
DataSets{1,58}='lux10_20141207T1540';
DataSets{1,59}='lux10_20141210T2154';
DataSets{1,60}='lux10_20141211T0125';
DataSets{1,61}='lux10_20141211T0501';
DataSets{1,62}='lux10_20141211T0834';
DataSets{1,63}='lux10_20141211T1202';
DataSets{1,64}='lux10_20141211T1501';
DataSets{1,65}='lux10_20150227T2215';
DataSets{1,66}='lux10_20150228T0549';
DataSets{1,67}='lux10_20150228T1323';
DataSets{1,68}='lux10_20150228T2057';
DataSets{1,69}='lux10_20150301T0151';
DataSets{1,70}='lux10_20150301T0930';
DataSets{1,71}='lux10_20150301T1710';
DataSets{1,72}='lux10_20150302T0047';
DataSets{1,73}='lux10_20150302T0827';
DataSets{1,74}='lux10_20150302T1600';
DataSets{1,75}='lux10_20150302T2335';
DataSets{1,76}='lux10_20150303T0707';
DataSets{1,77}='lux10_20150303T1439';
DataSets{1,78}='lux10_20150303T2216';
DataSets{1,79}='lux10_20150304T0549';
DataSets{1,80}='lux10_20150304T1325';
DataSets{1,81}='lux10_20150304T2059';
DataSets{1,82}='lux10_20150305T0441';
DataSets{1,83}='lux10_20150305T1214';
DataSets{1,84}='lux10_20150305T1946';
DataSets{1,85}='lux10_20150306T0324';
DataSets{1,86}='lux10_20150306T1343';
DataSets{1,87}='lux10_20150306T2122';
DataSets{1,88}='lux10_20150307T0504';
DataSets{1,89}='lux10_20150307T1342';
DataSets{1,90}='lux10_20150307T2043';
DataSets{1,91}='lux10_20150308T0326';
DataSets{1,92}='lux10_20150308T1108';
DataSets{1,93}='lux10_20150308T1854';
DataSets{1,94}='lux10_20150309T0231';
%}
DataSets{1,1402}='lux10_20160227T0653';
DataSets{1,1403}='lux10_20160227T1423';
DataSets{1,1404}='lux10_20160227T2155';
DataSets{1,1405}='lux10_20160228T0051';
DataSets{1,1406}='lux10_20160228T0822';
DataSets{1,1407}='lux10_20160228T1555';
DataSets{1,1408}='lux10_20160228T2335';
DataSets{1,1409}='lux10_20160229T0705';
DataSets{1,1410}='lux10_20160229T1438';
DataSets{1,1411}='lux10_20160229T2211';
DataSets{1,1412}='lux10_20160301T0545';
DataSets{1,1413}='lux10_20160301T1543';
DataSets{1,1414}='lux10_20160301T2322';
DataSets{1,1415}='lux10_20160302T0658';
DataSets{1,1416}='lux10_20160303T0052';
DataSets{1,1417}='lux10_20160303T0823';
DataSets{1,1418}='lux10_20160303T1600';
DataSets{1,1419}='lux10_20160303T2327';
DataSets{1,1420}='lux10_20160304T0654';
DataSets{1,1421}='lux10_20160304T1426';
DataSets{1,1422}='lux10_20160304T2156';
DataSets{1,1423}='lux10_20160305T0529';
DataSets{1,1424}='lux10_20160305T1303';
DataSets{1,1425}='lux10_20160305T2031';
DataSets{1,1426}='lux10_20160306T0051';
DataSets{1,1427}='lux10_20160306T0821';
DataSets{1,1428}='lux10_20160306T1547';
DataSets{1,1429}='lux10_20160306T2105';
DataSets{1,1430}='lux10_20160307T0434';
DataSets{1,1431}='lux10_20160307T1120';
DataSets{1,1432}='lux10_20160307T1857';
DataSets{1,1433}='lux10_20160308T0229';
DataSets{1,1434}='lux10_20160308T1004';
DataSets{1,1435}='lux10_20160308T1737';
DataSets{1,1436}='lux10_20160309T0024';
DataSets{1,1437}='lux10_20160309T0422';
DataSets{1,1438}='lux10_20160310T0003';
DataSets{1,1439}='lux10_20160310T0730';
DataSets{1,1440}='lux10_20160310T1504';
DataSets{1,1441}='lux10_20160310T2059';
DataSets{1,1442}='lux10_20160311T0436';
DataSets{1,1443}='lux10_20160311T1206';
DataSets{1,1444}='lux10_20160311T1939';
DataSets{1,1445}='lux10_20160312T0313';
DataSets{1,1446}='lux10_20160312T1047';
DataSets{1,1447}='lux10_20160312T1826';
DataSets{1,1448}='lux10_20160313T0051';
DataSets{1,1449}='lux10_20160313T0924';
DataSets{1,1450}='lux10_20160313T1659';
DataSets{1,1451}='lux10_20160314T0027';
DataSets{1,1452}='lux10_20160314T0803';
DataSets{1,1453}='lux10_20160314T1537';
DataSets{1,1454}='lux10_20160314T2310';
DataSets{1,1455}='lux10_20160315T0649';
DataSets{1,1456}='lux10_20160315T1059';
DataSets{1,1457}='lux10_20160315T1833';
DataSets{1,1458}='lux10_20160316T0207';
DataSets{1,1459}='lux10_20160316T0936';
DataSets{1,1460}='lux10_20160317T0902';
DataSets{1,1461}='lux10_20160317T1632';
DataSets{1,1462}='lux10_20160318T0002';
DataSets{1,1463}='lux10_20160318T0735';
DataSets{1,1464}='lux10_20160318T1502';
DataSets{1,1465}='lux10_20160318T2230';
DataSets{1,1466}='lux10_20160319T0554';
DataSets{1,1467}='lux10_20160319T1329';
DataSets{1,1468}='lux10_20160319T2059';
DataSets{1,1469}='lux10_20160320T0151';
DataSets{1,1470}='lux10_20160320T0918';
DataSets{1,1471}='lux10_20160320T1647';
DataSets{1,1472}='lux10_20160321T0020';
DataSets{1,1473}='lux10_20160321T0749';
DataSets{1,1474}='lux10_20160321T1518';
DataSets{1,1475}='lux10_20160321T2250';
DataSets{1,1476}='lux10_20160322T0622';
DataSets{1,1477}='lux10_20160322T1511';
DataSets{1,1478}='lux10_20160322T2245';
DataSets{1,1479}='lux10_20160323T0619';
DataSets{1,1480}='lux10_20160323T1344';
DataSets{1,1481}='lux10_20160323T2110';
DataSets{1,1482}='lux10_20160324T0440';
DataSets{1,1483}='lux10_20160324T2356';
DataSets{1,1484}='lux10_20160325T0722';
DataSets{1,1485}='lux10_20160325T1253';
DataSets{1,1486}='lux10_20160325T2026';
DataSets{1,1487}='lux10_20160326T0356';
DataSets{1,1488}='lux10_20160326T1126';
DataSets{1,1489}='lux10_20160326T1853';
DataSets{1,1490}='lux10_20160327T0151';
DataSets{1,1491}='lux10_20160327T0920';
DataSets{1,1492}='lux10_20160327T1649';
DataSets{1,1493}='lux10_20160328T0024';
DataSets{1,1494}='lux10_20160328T0751';
DataSets{1,1495}='lux10_20160328T1517';
DataSets{1,1496}='lux10_20160328T2246';
DataSets{1,1497}='lux10_20160329T0613';
DataSets{1,1498}='lux10_20160329T1529';
DataSets{1,1499}='lux10_20160329T2304';
DataSets{1,1500}='lux10_20160330T0635';
DataSets{1,1501}='lux10_20160331T0046';
DataSets{1,1502}='lux10_20160331T0816';
%{
DataSets{1,1503}='lux10_20160331T1547';
DataSets{1,1504}='lux10_20160331T2317';
DataSets{1,1505}='lux10_20160401T0650';
DataSets{1,1506}='lux10_20160401T1422';
DataSets{1,1507}='lux10_20160401T2151';
DataSets{1,1508}='lux10_20160402T0523';
DataSets{1,1509}='lux10_20160402T1257';
DataSets{1,1510}='lux10_20160402T2028';
DataSets{1,1511}='lux10_20160403T0151';
DataSets{1,1512}='lux10_20160403T0919';
DataSets{1,1513}='lux10_20160403T1652';
DataSets{1,1514}='lux10_20160404T0022';
DataSets{1,1515}='lux10_20160404T0750';
DataSets{1,1516}='lux10_20160404T1232';
DataSets{1,1517}='lux10_20160404T2000';
DataSets{1,1518}='lux10_20160405T0329';
DataSets{1,1519}='lux10_20160405T1102';
DataSets{1,1520}='lux10_20160405T1838';
DataSets{1,1521}='lux10_20160406T0213';
DataSets{1,1522}='lux10_20160406T1005';
DataSets{1,1523}='lux10_20160406T1739';
DataSets{1,1524}='lux10_20160407T0111';
DataSets{1,1525}='lux10_20160407T1110';
DataSets{1,1526}='lux10_20160407T2225';
DataSets{1,1527}='lux10_20160408T0551';
DataSets{1,1528}='lux10_20160408T1322';
DataSets{1,1529}='lux10_20160408T2238';
DataSets{1,1530}='lux10_20160409T0609';
DataSets{1,1531}='lux10_20160409T1344';
DataSets{1,1532}='lux10_20160409T2114';
DataSets{1,1533}='lux10_20160410T0151';
DataSets{1,1534}='lux10_20160410T0916';
DataSets{1,1535}='lux10_20160410T1650';
DataSets{1,1536}='lux10_20160411T0020';
DataSets{1,1537}='lux10_20160411T0757';
DataSets{1,1538}='lux10_20160411T1526';
DataSets{1,1539}='lux10_20160411T2256';
DataSets{1,1540}='lux10_20160412T0626';
DataSets{1,1541}='lux10_20160412T1329';
DataSets{1,1542}='lux10_20160412T2056';
DataSets{1,1543}='lux10_20160413T0428';
DataSets{1,1544}='lux10_20160413T0940';
DataSets{1,1545}='lux10_20160414T0424';
DataSets{1,1546}='lux10_20160414T1149';
DataSets{1,1547}='lux10_20160414T1855';
DataSets{1,1548}='lux10_20160415T0225';
DataSets{1,1549}='lux10_20160415T0914';
DataSets{1,1550}='lux10_20160415T1645';
DataSets{1,1551}='lux10_20160416T0022';
DataSets{1,1552}='lux10_20160416T0800';
DataSets{1,1553}='lux10_20160417T0151';
DataSets{1,1554}='lux10_20160417T0926';
DataSets{1,1555}='lux10_20160417T1659';
DataSets{1,1556}='lux10_20160418T0037';
DataSets{1,1557}='lux10_20160418T0808';
DataSets{1,1558}='lux10_20160418T1346';
DataSets{1,1559}='lux10_20160418T1718';
DataSets{1,1560}='lux10_20160419T0052';
DataSets{1,1561}='lux10_20160419T0942';
DataSets{1,1562}='lux10_20160419T1818';
DataSets{1,1563}='lux10_20160420T0156';
DataSets{1,1564}='lux10_20160420T0907';
DataSets{1,1565}='lux10_20160420T1642';
DataSets{1,1566}='lux10_20160421T0020';
DataSets{1,1567}='lux10_20160421T0750';
DataSets{1,1568}='lux10_20160421T1515';
DataSets{1,1569}='lux10_20160421T2247';
DataSets{1,1570}='lux10_20160422T0623';
DataSets{1,1571}='lux10_20160422T1346';
DataSets{1,1572}='lux10_20160422T1930';
DataSets{1,1573}='lux10_20160423T0305';
DataSets{1,1574}='lux10_20160423T1038';
DataSets{1,1575}='lux10_20160423T1805';
DataSets{1,1576}='lux10_20160424T0151';
DataSets{1,1577}='lux10_20160424T0927';
DataSets{1,1578}='lux10_20160424T1702';
DataSets{1,1579}='lux10_20160425T0035';
DataSets{1,1580}='lux10_20160425T0806';
DataSets{1,1581}='lux10_20160425T1541';
DataSets{1,1582}='lux10_20160425T2310';
DataSets{1,1583}='lux10_20160426T0637';
DataSets{1,1584}='lux10_20160426T1312';
DataSets{1,1585}='lux10_20160426T2036';
DataSets{1,1586}='lux10_20160427T0407';
DataSets{1,1587}='lux10_20160427T1136';
DataSets{1,1588}='lux10_20160428T0358';
DataSets{1,1589}='lux10_20160428T1120';
DataSets{1,1590}='lux10_20160428T1844';
DataSets{1,1591}='lux10_20160429T0212';
DataSets{1,1592}='lux10_20160429T0943';
DataSets{1,1593}='lux10_20160429T1708';
DataSets{1,1594}='lux10_20160430T0033';
DataSets{1,1595}='lux10_20160430T0805';
DataSets{1,1596}='lux10_20160430T1528';
DataSets{1,1597}='lux10_20160430T2255';
DataSets{1,1598}='lux10_20160501T0151';
DataSets{1,1599}='lux10_20160501T0922';
DataSets{1,1600}='lux10_20160501T1647';
DataSets{1,1601}='lux10_20160502T0012';
%}
addpath('/scratch/dkhaitan/IsolatedS1S2/r1');
RunSetting='/scratch/dkhaitan/IsolatedS1S2/r1/Settings_SR21_BH.xml';


cd '/scratch/dkhaitan/IsolatedS1S2/';
DataOutPath = ['/scratch/dkhaitan/IsolatedS1S2/output/'];
MaxDriftLength = 35000;


	DataSet=DataSets{1,iij};
	disp([DataSet])
	DataPath = ['/scratch/dkhaitan/LUXData/IsolatedS1S2/' DataSet '/'];
	FileList = dir([DataPath '*.dat']);
	FileListToRun = FileList;%datasample(FileList,NumberOfFiles,'Replace',false);

	Setting = XMLReader_framework(RunSetting);
	temp_isolated_s1_rate = 0;
	isolated_s2_rate = 0;
	s1_spect = zeros(1,1000);
	s2_spect = zeros(1,1000);
	livetime=0;
	counter=0;



	parfor ii=1:200
		try
			disp([ii,FileListToRun(ii).name]);
			[FileOut,ee_merge] = Level2(FileListToRun(ii).name,DataPath,DataOutPath,DataSet,Setting);
			[temp_isolated_s1_rate,temp_isolated_s2_rate,temp_s1_spect,temp_s2_spect,temp_livetime] = IsolateS1S2s(ee_merge,MaxDriftLength);
			%{
			isolated_s1_rate = isolated_s1_rate + temp_isolated_s1_rate;
			isolated_s2_rate = isolated_s1_rate + temp_isolated_s2_rate;
			s1_spect = s1_spect + temp_s1_spect
			s2_spect = s2_spect + temp_s2_spect
			livetime = livetime+temp_livetime
			%}
			%{
			isolated_s1_rate = (livetime_total*isolated_s1_rate + livetime*temp_isolated_s1_rate)/(livetime_total + livetime);
			isolated_s2_rate = (livetime_total*isolated_s2_rate + livetime*temp_isolated_s2_rate)/(livetime_total + livetime);
			s1_spect = (s1_spect*livetime_total + livetime*temp_s1_spect)/(livetime_total + livetime);
			s2_spect = (s2_spect*livetime_total + livetime*temp_s2_spect)/(livetime_total + livetime);
			livetime_total=livetime_total + livetime;
			%}
			%disp([ii,temp_isolated_s1_rate,temp_isolated_s2_rate,temp_livetime]);
			fid_s2_spect=fopen([DataOutPath DataSet 's2_spect.csv'],'at');
			fprintf(fid_s2_spect,'%s,',num2str(temp_isolated_s2_rate));
			fprintf(fid_s2_spect,'%s,',num2str(temp_livetime));
			dlmwrite([DataOutPath DataSet 's2_spect.csv'],temp_s2_spect,'-append','delimiter',',');
			fclose(fid_s2_spect);

			fid_s1_spect=fopen([DataOutPath DataSet 's1_spect.csv'],'at');
			fprintf(fid_s1_spect,'%s,',num2str(temp_isolated_s1_rate));
			fprintf(fid_s2_spect,'%s,',num2str(temp_livetime));
			dlmwrite([DataOutPath DataSet 's1_spect.csv'],temp_s1_spect,'-append','delimiter',',');
			fclose(fid_s1_spect);
		catch exception
			disp(['Fail ' FileList(ii).name])
		end
		%clearvars -global -except DataSets DataSet DataPath DataOutPath MaxDriftLength FileList FileListToRun NumberOfFiles Setting isolated_s1_rate isolated_s2_rate s1_spect s2_spect livetime counter
	end
	dummy = csvread([DataOutPath DataSet 's2_spect.csv']);
	dummy = sum(dummy);
	isolated_s2_rate = dummy(1)/dummy(2);
	s2_spect = dummy(3:length(dummy))/dummy(2);
	dummy = csvread([DataOutPath DataSet 's1_spect.csv']);
	dummy = sum(dummy);
	isolated_s1_rate = dummy(1)/dummy(2);
	s1_spect = dummy(3:length(dummy))/dummy(2);

	fid_s2_spect=fopen([DataOutPath 'total_s2_spect.csv'],'at');
	fprintf(fid_s2_spect,'%s,',DataSet);
	fprintf(fid_s2_spect,'%s,',num2str(isolated_s2_rate));
	fprintf(fid_s2_spect,'%s,',num2str(dummy(2)));
	dlmwrite([DataOutPath 'total_s2_spect_real1.csv'],s2_spect,'-append','delimiter',',');
	fclose(fid_s2_spect);

	fid_s1_spect=fopen([DataOutPath 'total_s1_spect.csv'],'at');
	fprintf(fid_s1_spect,'%s,',DataSet);
	fprintf(fid_s1_spect,'%s,',num2str(isolated_s1_rate));
	fprintf(fid_s2_spect,'%s,',num2str(dummy(2)));
	dlmwrite([DataOutPath 'total_s1_spect_real1.csv'],s1_spect,'-append','delimiter',',');
	fclose(fid_s1_spect);
	delete([DataOutPath DataSet 's2_spect.csv'])
	delete([DataOutPath DataSet 's1_spect.csv'])
	
end
