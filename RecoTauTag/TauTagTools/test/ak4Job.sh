mkdir -p ${dir}/job_${3}
cd ${dir}/job_${3}
currDir=`pwd`
echo 'Now working in '${currDir}
eval `scramv1 runtime -sh`
cp ../runDiff_cfg.py .
cmsRun runDiff_cfg.py maxEvents=${1} skipEvents=${2}
