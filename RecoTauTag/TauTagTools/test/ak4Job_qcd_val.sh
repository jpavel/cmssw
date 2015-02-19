cd ${dir}
currDir=`pwd`
echo 'Now working in '${currDir}
eval `scramv1 runtime -sh`
cmsRun runDiff_qcd_val_cfg.py maxEvents=-1
