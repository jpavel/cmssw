import os,sys
argList=sys.argv
if len(argList) != 4:
  print 'Wrong number of arguments!'
  sys.exit(1)

jobName=argList[1]
eventsPerJob=argList[2]
numberOfJobs=argList[3]
batchScript='ak4Job.sh'
cfgScript='runDiff_cfg.py'
currDir=os.getcwd()

print 'This job will run from dir '+jobName
print 'The number of events per job will be '+str(eventsPerJob)
print 'The number of jobs will be '+str(numberOfJobs)
cont=raw_input('Do you agree? (Y/N) ')

if cont!='Y' :
  print 'Ending execution!'
  sys.exit(1)

from subprocess import call
call(["mkdir","-p",jobName]) # make dir if not exist
#call(["cp",cfgScript,jobName])

f=open(batchScript)
batchCommands=f.readlines()
f.close()

# cmd_f=open(jobName+'/'+batchScript,'w')
# cmd_f.write('dir='+currDir+'/'+jobName+'\n')
# for line in batchCommands:
#   cmd_f.write(line)

# cmd_f.close()
# skipEvents=0
for j_number in range(int(numberOfJobs)):
  print j_number
  skipEvents=j_number*int(eventsPerJob)
  call(["mkdir","-p",jobName+"/job_"+str(j_number)])
  command=currDir+'/'+jobName+'/job_'+str(j_number)+'/'+batchScript+' '+str(eventsPerJob)+' '+str(skipEvents)+' '+str(j_number)
  print command
  call(["cp",cfgScript,jobName+"/job_"+str(j_number)])
  os.chdir(jobName+"/job_"+str(j_number))
  print os.getcwd()
  cmd_f=open(batchScript,'w')
  cmd_f.write('dir='+currDir+'/'+jobName+'/job_'+str(j_number)+'\n')
  for line in batchCommands:
    cmd_f.write(line)
  cmd_f.close()
  call(["chmod","+x",batchScript])
  call(["bsub","-q","1nh",command])
  os.chdir(currDir)

print os.getcwd()
  #  bsub -q test -W 1:00 '/afs/cern.ch/work/j/jez/private/CMS/validation/140402_AK4FRproblem/ak5/CMSSW_7_1_0_pre4/src/test_batch/ak4Job.sh 1000 0 1
  #submission
  #call(["bsub","-q","8nm","-W","0:01",command])
#  call(["bsub","-q","8nm",command])
