import os,sys
argList=sys.argv
if len(argList) != 4:
  print 'Wrong number of arguments! Syntax is python prepareJob_splitByFiles.py <jobName> <filesPerJob> <numberOfJobs>'
  sys.exit(1)

jobName=argList[1]
filesPerJob=argList[2]
numberOfJobs=argList[3]
batchScript='ak4Job_ztt.sh'
cfgScript='runDiff_ztt_val_cfg.py'
cfgScript_part1='runDiff_ztt_val_cfg_part1.py'
cfgScript_part2='runDiff_ztt_val_cfg_part2.py'
filelist='file.list_fullsim'
currDir=os.getcwd()


print 'This job will run from dir '+jobName
print 'The number of files per job will be '+str(filesPerJob)
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

g=open(filelist)
listOfFiles=g.readlines()
g.close()

h1=open(cfgScript_part1)
cfgCommands1=h1.readlines()
h1.close()

h2=open(cfgScript_part2)
cfgCommands2=h2.readlines()
h2.close()


# cmd_f=open(jobName+'/'+batchScript,'w')
# cmd_f.write('dir='+currDir+'/'+jobName+'\n')
# for line in batchCommands:
#   cmd_f.write(line)

# cmd_f.close()
# skipEvents=0
for j_number in range(int(numberOfJobs)):
  print j_number
  firstFile=j_number*int(filesPerJob)
  lastFile=firstFile+int(filesPerJob)
  if lastFile > len(listOfFiles):
    lastFile=len(listOfFiles)
  print 'files '+str(firstFile)+' to '+str(lastFile)
  call(["mkdir","-p",jobName+"/job_"+str(j_number)])
  command=currDir+'/'+jobName+'/job_'+str(j_number)+'/'+batchScript
  print command
  os.chdir(jobName+"/job_"+str(j_number))
  print os.getcwd()
  cmd_f=open(batchScript,'w')
  cmd_f.write('dir='+currDir+'/'+jobName+'/job_'+str(j_number)+'\n')
  for line in batchCommands:
    cmd_f.write(line)
  cmd_f.close()
  cfg_f=open(cfgScript,'w')
  for line in cfgCommands1:
    cfg_f.write(line)
  for lineNum in range(firstFile,lastFile):
    cfg_f.write(listOfFiles[lineNum])
  for line in cfgCommands2:
    cfg_f.write(line)
  cfg_f.close()
  call(["chmod","+x",batchScript])
  call(["bsub","-q","1nh",command])
  os.chdir(currDir)

print os.getcwd()
  #  bsub -q test -W 1:00 '/afs/cern.ch/work/j/jez/private/CMS/validation/140402_AK4FRproblem/ak5/CMSSW_7_1_0_pre4/src/test_batch/ak4Job.sh 1000 0 1
  #submission
  #call(["bsub","-q","8nm","-W","0:01",command])
#  call(["bsub","-q","8nm",command])
