# Execute all NEW H2 identification/robustness probes; keep individual output files.
# Place this file, the six JPR_H2_CAUSAL_*.R files, the earlier
# JPR_H2_LWJ_00_COMMON.R and two LWJ CSVs in Windows Downloads.
BASE<-'C:/Users/USERNAME/Downloads'
SCRIPTS<-c('JPR_H2_CAUSAL_07_SAMPLE_AUDIT.R',
 'JPR_H2_CAUSAL_08_WITHIN_DISTRICT_MONTH.R',
 'JPR_H2_CAUSAL_09_FROZEN_2017_EVENT_STUDY.R',
 'JPR_H2_CAUSAL_10_FROZEN_SNAPSHOT_SENSITIVITY.R',
 'JPR_H2_CAUSAL_11_CONSOLIDATION_AND_CLUSTERING.R',
 'JPR_H2_CAUSAL_12_JOINT_ACLED_LWJ.R')
LOG<-data.frame(script=character(),status=character(),message=character())
for(s in SCRIPTS){
 cat('\n==== RUNNING',s,'====\n')
 p<-file.path(BASE,s)
 outcome<-tryCatch({source(p,local=new.env(parent=globalenv()));list(status='ok',message='')},
                   error=function(e)list(status='error',message=conditionMessage(e)))
 LOG<-rbind(LOG,data.frame(script=s,status=outcome$status,message=outcome$message))
 cat('RESULT:',outcome$status,outcome$message,'\n')
}
print(LOG,row.names=FALSE)
write.csv(LOG,file.path(BASE,'JPR_H2_CAUSAL_RUN_LOG.csv'),row.names=FALSE)
cat('Diagnostics are not a substitute for an exogenous surveillance shock.\n')
