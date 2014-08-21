cd ${dir}
currDir=`pwd`
echo 'Now working in '${currDir}
eval `scramv1 runtime -sh`
cmsRun runDiff_qcd_cfg.py maxEvents=${1} skipEvents=${2}
