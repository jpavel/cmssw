1) Upload the payloads to the PREP dropbox:

python uploadConstants.py --version=2014Apr29 --sqlite_file=RecoTauTag_MVAs_2014Apr29.db --prep 


2) Now the files can be checked here:

https://cms-conddb-prod.cern.ch/logs/dropBox/


3) Check the uploaded test payloads. 

cmscond_list_iov -c frontier://FrontierPrep/CMS_COND_PAT_001 -a | grep RecoTauTag


NOTE: If that doesn't work, try 
cmscond_list_iov -P /afs/cern.ch/cms/DB/conddb -c oracle://cms_orcoff_prep/CMS_COND_PAT_001 -a | grep RecoTauTag



4) Now run the "uploadConstants.py" script again. Ommit the "--prep" option this time

python uploadConstants.py --version=2014Apr29 --sqlite_file=RecoTauTag_MVAs_2014Apr29.db
 

