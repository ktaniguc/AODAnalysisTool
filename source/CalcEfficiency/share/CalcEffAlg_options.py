if 'LocalInputFileList' in locals():
    print "LocalInputFileList is already set"
else:
    #LocalInputFileList="./list/Wplusmunu_isMoreCandtest.list"
    #LocalInputFileList="./list/JpsiCollimated_isMoreCandtest.list"
    #LocalInputFileList="./list/testNSW.list"
    #LocalInputFileList="./list/Jpsimu4mu20_all.list"
    #LocalInputFileList="./list/Zmumu_350files.list"
    #LocalInputFileList="test_ynoguchi.list"
    LocalInputFileList="./list/trigArt_0402.list"
print "LocalInputFileList is"
print LocalInputFileList
f = open( LocalInputFileList, 'r' )
InputFileList = f.read().splitlines()
print InputFileList
#-----------------------------------------------------------------------------
# Athena imports
#-----------------------------------------------------------------------------
from AthenaCommon.Constants import *
from AthenaCommon.AppMgr import theApp
from AthenaCommon.AppMgr import ServiceMgr
from AthenaCommon.AppMgr import ToolSvc
#import AthenaPoolCnvSvc.ReadAthenaPool
from AthenaPoolCnvSvc.PoolReadConfig import PoolReadCfg
from AthenaCommon.AlgSequence import AlgSequence
job = AlgSequence()
# Setup the Run III behavior    

#-----------------------------------------------------------------------------
# Message Service
#-----------------------------------------------------------------------------
# Set output level threshold (2=DEBUG, 3=INFO, 4=WARNING, 5=ERROR, 6=FATAL )
ServiceMgr.MessageSvc.OutputLevel = ERROR
import AthenaServices
AthenaServices.AthenaServicesConf.AthenaEventLoopMgr.OutputLevel = ERROR

#-----------------------------------------------------------------------------
# Input Datasets
#-----------------------------------------------------------------------------
from AthenaCommon.AthenaCommonFlags import athenaCommonFlags as acf
if not acf.EvtMax.is_locked():
  acf.EvtMax=-1
acf.FilesInput = InputFileList
#acf.FilesInput += [
#"/gpfs/fs6001/yyazawa/data/valid1.424100.Pythia8B_A14_CTEQ6L1_Jpsimu4mu4.recon.AOD.e5112_s3091_r9122_tid10750758_00/AOD.10750758._000194.pool.root.1",
#"/gpfs/fs6001/yyazawa/data/valid1.424100.Pythia8B_A14_CTEQ6L1_Jpsimu4mu4.recon.AOD.e5112_s2887_r9026_tid10522817_00/AOD.10522817._000212.pool.root.1"
#]
#-----------------------------------------------------------------------------
# Algorithms
#-----------------------------------------------------------------------------
rec.doCBNT=False
from RecExConfig.RecFlags import rec
rec.doTrigger=True
from RecExConfig.RecAlgsFlags import recAlgs
recAlgs.doTrigger=True
recAlgs.doAtlfast=False
recAlgs.doMonteCarloReact=False
from TriggerJobOpts.TriggerFlags import TriggerFlags
#TriggerFlags.doTriggerConfigOnly=True
print("TriggerFlags.doMT = true")
TriggerFlags.doMT=True
print("TriggerFlags.doMT = true end")

rec.doWriteAOD=False
rec.doWriteESD=False
rec.doWriteTAG=False
rec.doAOD=False
rec.doDPD=False
rec.doESD=False
doTAG=False
rec.doTruth=False
rec.doRecoTiming=False
#rec.doDetStatus=False
rec.doShowSizeStatistics=False
rec.readTAG=False
rec.readRDO=False
rec.doHist=False
rec.doContainerRemapping=False
rec.doJiveXML=False
rec.doEdmMonitor=False
rec.doDumpPoolInputContent=False
rec.doHeavyIon=False
rec.doHIP=False
rec.doWriteBS=False
rec.doPhysValMonHists=False
rec.doVP1=False
rec.doJiveXML=False
rec.doMuon=False
rec.doCheckDictionary=False
rec.doFileMetaData=False
rec.doCalo=False
rec.doAODCaloCells=False
rec.doEgamma=False

## Set the Athena configuration flags
#from AthenaConfiguration.AllConfigFlags import ConfigFlags
#from AthenaConfiguration.AutoConfigFlags import GetFileMD
#
#from AthenaMonitoring.DQConfigFlags import allSteeringFlagsOff
#allSteeringFlagsOff()
## bring things into scope
#from AthenaMonitoring.DQConfigFlags import allSteeringFlagsOff
#ConfigFlags.dump()
#
#ConfigFlags.lock()
#
## Initialize configuration object, add accumulator, merge, and run.
#from AthenaConfiguration.MainServicesConfig import MainServicesCfg
## load DQ
#
#rec.doESD.set_Value_and_Lock(False) # uncomment if do not run ESD making algorithms
#rec.doWriteESD.set_Value_and_Lock(False) # uncomment if do not write ESD
#rec.doAOD.set_Value_and_Lock(False) # uncomment if do not run AOD making algorithms
#rec.doWriteAOD.set_Value_and_Lock(False) # uncomment if do not write AOD
#rec.doWriteTAG.set_Value_and_Lock(False) # uncomment if do not write TAG
#
#include("RecExCommon/RecExCommon_topOptions.py")
#
########
#ToolSvc.TrigDecisionTool.TrigDecisionKey='xTrigDecision'
#from TrigDecisionTool.TrigDecisionToolConf import Trig__TrigDecisionTool
#from TrigNavigation.TrigNavigationConf import HLT__Navigation
#ToolSvc += CfgMgr.Trig__TrigDecisionTool( "TrigDecisionTool" )
#ToolSvc += CfgMgr.HLT__Navigation( "Navigation" )
#ToolSvc.TrigDecisionTool.Navigation.ReadonlyHolders=True
########
#ServiceMgr.MessageSvc.setWarning += [ "", "HolderFactory" ]
#ServiceMgr.MessageSvc.infoLimit = 99999999
#ServiceMgr.MessageSvc.debugLimit = 99999999

# GRL
ToolSvc += CfgMgr.GoodRunsListSelectionTool("MyGRLTool",GoodRunsListVec=["current_grl.xml"])

include("RecExCommon/RecExCommon_topOptions.py")

from CalcEfficiency.CalcEfficiencyConf import *
#
#preparator = CalcEffAlg__TagAndProbe()
#preparator.OutputLevel = DEBUG
#ToolSvc += preparator

#ServiceMgr.CalcEffAlg.OutputLevel = INFO
#ServiceMgr.MessageSvc.debugLimit = 99999999
#ServiceMgr.MessageSvc.defaultLimit = 99999999

#/MufastHypoConfig.OutputLevel = DEBUG
#MucombHypoConfig.OutputLevel = DEBUG
#AthToolSupport/AsgTools
#from AthToolSupport import AsgTools
#AsgTools.OutputLevel = DEBUG
#AsgTools_UserCode.OutputLevel = DEBUG

#from TrigDecisionTool.TrigDecisionToolConf import Trig__TrigDecisionTool
#job+=Trig__TrigDecisionTool
#job.Trig__TrigDecisionTool.OutputLevel = DEBUG

#from TrigDecisionMaker.TrigDecisionMakerConfig import TriggerBitsMakerTool
#job+=TriggerBitsMakerTool()
#job.TriggerBitsMakerTool.OutputLevel = DEBUG
#####
#from TrigDecisionMaker.TrigDecisionMakerConfig import TrigDecisionMakerMT
#job+=TrigDecisionMakerMT()
#job.TrigDecMakerMT.OutputLevel = DEBUG

# Setup logs
#from AthenaCommon.Logging import log
#from AthenaCommon.Constants import *
#log.setLevel(DEBUG)


job += CalcEffAlg( message = 1,
                   OutputLevel = DEBUG,
                   OutputFile = "/gpfs/fs7001/ktaniguc/outputfile/OutputCalcEff/test0621.root",
                   #TapMethod = "NoMass",
                   #TapMethod = "Jpsi",
                   TapMethod = "Ztap",
                   Extrapolate = True,
                   GRL = False,
                   DataType = "data17"
                 )
                   #GRL = True
                   #OutputLevel = ERROR,
                   #TapMethod = "NoTag",
                   #TapMethod = "Ztap",
                   #TapMethod = "JPZtap",
                   #OutputFile = "test_inEITrkConfFalsev3.root",
                   #Datatype = "data16"

include("TriggerTest/TriggerTestCommon.py")
                     
print job

#-----------------------------------------------------------------------------
