#! /bin/sh
DATE=$(date '+%Y_%m_%d_%H_%M')

mkdir run_$DATE

echo "#! /bin/sh" > run_$DATE/run.sh
echo "" >> run_$DATE/run.sh
echo "export ATLAS_LOCAL_ROOT_BASE=/cvmfs/atlas.cern.ch/repo/ATLASLocalRootBase" >> run_$DATE/run.sh
echo "alias setupATLAS='source ${ATLAS_LOCAL_ROOT_BASE}/user/atlasLocalSetup.sh'" >> run_$DATE/run.sh
echo "setupATLAS" >> run_$DATE/run.sh
echo "" >> run_$DATE/run.sh
echo "cd $TestArea/../" >> run_$DATE/run.sh
echo "source ./setup.sh" >> run_$DATE/run.sh
echo "source $TestArea/*/setup.sh" >> run_$DATE/run.sh
echo 'export FRONTIER_SERVER="(serverurl=http://atlasfrontier-ai.cern.ch:8000/atlr)(serverurl=http://atlasfrontier2-ai.cern.ch:8000/atlr)(serverurl=http://atlasfrontier1-ai.cern.ch:8000/atlr)(serverurl=http://frontier-atlas.lcg.triumf.ca:3128/ATLAS_frontier)(serverurl=http://frontier-atlas1.lcg.triumf.ca:3128/ATLAS_frontier)(serverurl=http://frontier-atlas2.lcg.triumf.ca:3128/ATLAS_frontier)(serverurl=http://ccfrontier.in2p3.fr:23128/ccin2p3-AtlasFrontier)(serverurl=http://ccfrontier01.in2p3.fr:23128/ccin2p3-AtlasFrontier)(proxyurl=http://conddb-px01.icepp.jp:3128)(proxyurl=http://conddb-px02.icepp.jp:3128)(proxyurl=http://atlasbpfrontier.cern.ch:3127)(proxyurl=http://atlasbpfrontier.fnal.gov:3127)"' >> run_$DATE/run.sh
echo cd $TestArea/../run_bsub/run_$DATE/ >> run_$DATE/run.sh
echo "cp ../CalcEffAlg_optionsMT.py ./" >> run_$DATE/run.sh
echo "cp ../current_grl.xml ./" >> run_$DATE/run.sh
echo 'athena.py CalcEffAlg_optionsMT.py' >> run_$DATE/run.sh


chmod 744 run_$DATE/run.sh

cd run_$DATE
bsub -q 12h -n 2 -o ../LOGS/log_${DATE}.txt ./run.sh
cd ..
