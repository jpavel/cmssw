* exporting to PREP
  -> it is exported to oracle, in about 10-15 mins it appears in frontier

cmscond_export_iov -P /afs/cern.ch/cms/DB/conddb -s sqlite_file:RecoTauTag_MVAs_2014Jul07.db -d oracle://cms_orcoff_prep/CMS_COND_PAT_001 -i RecoTauTag_againstMuonMVAv1 -t RecoTauTag_againstMuonMVAv1

* to check

cmscond_list_iov -c frontier://FrontierPrep/CMS_COND_PAT_001 -t RecoTauTag_againstMuonMVAv1

or
cmscond_list_iov -c oracle://cms_orcoff_prep/CMS_COND_PAT_001 -P /afs/cern.ch/cms/DB/conddb -t RecoTauTag_againstMuonMVAv1

