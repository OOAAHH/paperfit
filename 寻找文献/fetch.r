library(easyPubMed)
#Use R env R432
#key1
organs <- c('Blood', 'Nose', 'Mouth', 'Pharynx', 'Larynx', 'Trachea', 'Bronchi', 'Lungs', 'Lung', 'Alveoli', 'Diaphragm', 'Respiratory system', 'Asthma', 'Bronchitis', 'Chronic Obstructive Pulmonary Disease', 'Pneumonia', 'Tuberculosis', 'Lung Cancer', 'Pulmonary Fibrosis', 'Pulmonary Hypertension', 'Cystic Fibrosis', 'Acute Respiratory Distress Syndrome', 'Pleural Effusion', 'Pneumothorax', 'Sarcoidosis', 'Lung Abscess', 'Flu', 'Influenza', 'COVID-19', 'SARS-cov-2', 'SARS', 'Allergy', 'Asthma', 'Silicosis', 'Pneumoconiosis', 'Pulmonary Fibrosis', 'COPD', 'Chronic Obstructive Pulmonary Disease', 'ILD', 'Interstitial Lung Disease', 'IPF', 'Idiopathic Pulmonary Fibrosis', 'KD', 'Kawasaki Disease', 'LSCC', 'Lung Squamous Cell Carcinoma', 'LUAD', 'Lung Adenocarcinoma', 'NSCLC', 'Non-Small Cell Lung Cancer', 'NSIP', 'Nonspecific Interstitial Pneumonia', 'SCC', 'Squamous Cell Carcinomas', 'SCLC', 'Small Cell Lung Cancer', 'SSc-ILD', 'Systemic Sclerosis-Associated Interstitial Lung Disease', 'cHP', 'Chronic Hypersensitivity Pneumonitis')
#key2
data_t=c('scRNA-seq','scATAC-seq','snRNA-seq','singlecell RNA sequencing',
'spatial','stero-seq','visium',
'BTR','TCR'
)
all_query=list()
for(i in organs){
all_query[[i]]=paste("(",i,")"," ","AND"," ","(",data_t,")",sep="")}
all_query1=unlist(all_query)###转变为向量
new_query=get_pubmed_ids(all_query1[1])
fetched_data=fetch_pubmed_data(new_query, encoding = "ASCII")
organs_PM_df=table_articles_byAuth(pubmed_data = fetched_data, 
included_authors = "first", 
max_chars = 0, 
encoding = "ASCII")
organs_PM_df['name']=paste(organs_PM_df$lastname,'et al',organs_PM_df$jabbrv)
organs_PM_df=organs_PM_df[,c('doi','name','year','title')]
colnames(organs_PM_df)=c('DOI', 'Name', 'Year', 'Title')

for(i in 2:length(all_query1)){
  new_query=get_pubmed_ids(all_query1[i])
  if(new_query$Count!=0){
  fetched_data=fetch_pubmed_data(new_query, encoding = "ASCII")
  new_PM_df=table_articles_byAuth(pubmed_data = fetched_data, 
  included_authors = "first", 
  max_chars = 0, 
  encoding = "ASCII")
  new_PM_df['name']=paste(new_PM_df$lastname,'et al',new_PM_df$jabbrv)
  new_PM_df=new_PM_df[,c('doi','name','year','title')]
  colnames(new_PM_df)<-c('DOI', 'Name', 'Year', 'Title')
  organs_PM_df=rbind(organs_PM_df,new_PM_df)}
}
# 保存organs_PM_df为CSV文件
write.csv(organs_PM_df, "~/filename.csv", row.names = FALSE)
