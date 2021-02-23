# Declare the package name:
atlas_subdir( CalcEfficiency )

# External dependencies:
find_package( CLHEP )
find_package( ROOT COMPONENTS GaudiInterface Core Tree MathCore Hist RIO pthread )

atlas_add_component( CalcEfficiency
#  src/*.cxx
  src/components/*.cxx
  INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS}
  LINK_LIBRARIES CalcEfficiencyLib )

atlas_add_library( CalcEfficiencyLib
  src/*.cxx
  PUBLIC_HEADERS CalcEfficiency 
  INCLUDE_DIRS ${ROOT_INCLUDE_DIRS} ${CLHEP_INCLUDE_DIRS}
  LINK_LIBRARIES ${ROOT_LIBRARIES} ${CLHEP_LIBRARIES} AthenaBaseComps xAODTrigMuon xAODTrigger StoreGateLib EventInfo GaudiKernel xAODEventInfo xAODMuon TrkParameters TrkVKalVrtFitterLib AsgTools GoodRunsListsLib TrigDecisionToolLib TrigConfInterfaces TrigConfxAODLib )

# Install files from the package:
atlas_install_headers( CalcEfficiency )
atlas_install_joboptions( run/CalcEffAlg_options.py )
