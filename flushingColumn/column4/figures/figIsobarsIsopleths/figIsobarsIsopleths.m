% extract and plot isobars and isopleths on a dissolved CO2 vs dissolved water diagram from solwcad output file; plot column evolution towards equilibrium as flushing proceeds

% Chiara P. Montagna, INGV Pisa, 2022/3
% chiara.montagna@ingv.it

% ============== INPUT ===============================
% solwcad output file 'file.out' contains columns:
% 1) pressure
% 3) total water
% 4) carbon_dioxide
% 5) dissolved_h2o_onliq 
% 6) dissolved_co2_onliq
% 8) exsolved CO2 on gas, wt frac

resultsFile = 'solwcadResultsIsopleth.out';

% pressures for isobars, MPa
pressures = [100 200 300 400];

% CO2 on gas for isopleths, wt frac
CO2s = [0.85 0.97];

% flushing evolution dir
flushingPath = '/disk/science/my_science/magmaDynamics/chamber/CO2flux/stromboli/solwcadColumn/pyFlushing/flushingColumn/column4/pyResults/';
% flushing steps to be plotted
inputSteps = [0 0.10 0.50 1.00 3.00 5.00 10.00];

% =============================================

% load data
data = load(resultsFile);

p = data(:,1);
totalWater = data(:,3);
totalCO2 = data(:,4);
dissWater = data(:,5);
dissCO2 = data(:,6);
exCO2 = data(:,8);

figure(1)
hold on

% loop on pressures
for i = 1:length(pressures)
    pres = pressures(i)*1.e6;
    indexP = find(p == pres);
    dissWaterP = dissWater(indexP);
    dissCO2P = dissCO2(indexP);
    dissCO2Sat = find(totalCO2(indexP) - dissCO2(indexP) > 0);
    dissWaterC = dissWaterP(dissCO2Sat);
    dissCO2C = dissCO2P(dissCO2Sat);

    figure(1)
    plot(dissWaterC, dissCO2C,'+')
end

for i = 1:length(CO2s)
    indexC = find(abs(exCO2 - CO2s(i)) < 0.001);
    plot(dissWater(indexC),dissCO2(indexC),'*')
end

ylim([0 7e-3])         
xlim([0 0.04])         
xlabel('Dissolved H2O')
ylabel('Dissolved CO2')

% plot flushing evolution
if inputSteps(1) == 0
    init = importdata([flushingPath 'initialColumnPy.out'],' ',2);
    plot(init.data(:,5),init.data(:,6),'lineWidth',2)
    start = 2;
else
    start = 1;
end

for i = start:length(inputSteps)
    in01 = importdata([flushingPath 'pythonColumn' num2str(inputSteps(i),'%.2f') '.out'],' ',2);
    plot(in01.data(:,5),in01.data(:,6),'lineWidth',2);
end