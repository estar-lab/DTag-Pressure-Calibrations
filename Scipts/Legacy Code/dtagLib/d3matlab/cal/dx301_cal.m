%Cal file for dx301, Dtag3X Core version
DEV.ID='c123251c';
DEV.NAME='DX301';
DEV.BUILT=[2016 4 28];
DEV.BUILDER='SF';
DEV.HAS={'stereo audio'};
BBFILE = ['badblocks_' DEV.ID(1:4) '_' DEV.ID(5:8) '.txt'] ;
try,
   DEV.BADBLOCKS = readbadblocks(['c:/mystuff/d3/host/d3usb/' BBFILE]) ;
catch,
   fprintf(' No bad block file\n') ;
end

TEMPR = struct ;
TEMPR.TYPE='ntc thermistor';
TEMPR.USE='conv_ntc';
TEMPR.UNIT='degrees Celcius';
TEMPR.METHOD='none';

BATT = struct ;
BATT.POLY=[6 0] ;
BATT.UNIT='Volt';

PRESS=struct;
PRESS.POLY=[31.628 2771.632 -77.709];
PRESS.METHOD='rough';
PRESS.LASTCAL=[2016 4 28];
PRESS.TREF = 20 ;
PRESS.UNIT='meters H20 salt';
PRESS.TC.POLY=[0];
PRESS.TC.SRC='BRIDGE';
PRESS.BRIDGE.NEG.POLY=[3 0];
PRESS.BRIDGE.NEG.UNIT='Volt';
PRESS.BRIDGE.POS.POLY=[6 0];
PRESS.BRIDGE.POS.UNIT='Volt';
PRESS.BRIDGE.RSENSE=200;
PRESS.BRIDGE.TEMPR.POLY=[314.0 -634.7] ;
PRESS.BRIDGE.TEMPR.UNIT='degrees Celcius';

ACC=struct;
ACC.TYPE='MEMS accelerometer';
ACC.POLY=[10.065 -4.959; 10.056 -4.983; 10.203 -5.065];
ACC.UNIT='g';
ACC.TREF = 20 ;
ACC.TC.POLY=[0; 0; 0];
ACC.PC.POLY=[0; 0; 0];
ACC.XC=zeros(3);
ACC.MAP=[-1 0 0;0 1 0;0 0 1];
ACC.MAPRULE='front-right-down';
ACC.METHOD='flips';
ACC.LASTCAL=[2016 4 28];

MAG=struct;
MAG.TYPE='magnetoresistive bridge';
MAG.POLY=[757.808 -257.04; 739.679 -224.252; 835.637 -159.642];
MAG.UNIT='Tesla';
MAG.TREF = 20 ;
MAG.TC.POLY=[0;0;0];
MAG.TC.SRC='BRIDGE';
MAG.PC.POLY=[0;0;0];
MAG.XC=zeros(3);
MAG.MAP=[0 1 0;1 0 0;0 0 1];
MAG.MAPRULE='front-right-down';
MAG.METHOD='';
MAG.LASTCAL=[2016 4 28];
MAG.BRIDGE.NEG.POLY=[3 0];
MAG.BRIDGE.NEG.UNIT='Volt';
MAG.BRIDGE.POS.POLY=[6 0];
MAG.BRIDGE.POS.UNIT='Volt';
MAG.BRIDGE.RSENSE=20;
MAG.BRIDGE.TEMPR.POLY=[541.91 -459.24] ;
MAG.BRIDGE.TEMPR.UNIT='degrees Celcius';

CAL=struct ;
CAL.TEMPR=TEMPR;
CAL.BATT=BATT;
CAL.PRESS=PRESS;
CAL.ACC=ACC;
CAL.MAG=MAG;

DEV.CAL = CAL ;
writematxml(DEV,'DEV','dx301.xml')

