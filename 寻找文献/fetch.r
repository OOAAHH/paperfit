#devtools::install_github("dami82/easyPubMed")
library(easyPubMed)
#Use R env R432
#key1
organs <- c('Blood','Nose','Mouth','Pharynx','Larynx','Trachea','Bronchi','Lungs','Lung','Alveoli','Diaphragm','Respiratory system','Asthma','Bronchitis','Chronic Obstructive Pulmonary Disease','Pneumonia','Tuberculosis','Lung Cancer','Pulmonary Fibrosis','Pulmonary Hypertension','Cystic Fibrosis','Acute Respiratory Distress Syndrome','Pleural Effusion','Pneumothorax','Sarcoidosis','Lung Abscess','Flu','Influenza','COVID-19','SARS-cov-2','SARS','Allergy','Asthma','Silicosis','Pneumoconiosis','Pulmonary Fibrosis','COPD','Chronic Obstructive Pulmonary Disease','ILD','Interstitial Lung Disease','IPF','Idiopathic Pulmonary Fibrosis','KD','Kawasaki Disease','LSCC','Lung Squamous Cell Carcinoma','LUAD','Lung Adenocarcinoma','NSCLC','Non-Small Cell Lung Cancer','NSIP','Nonspecific Interstitial Pneumonia','SCC','Squamous Cell Carcinomas','SCLC','Small Cell Lung Cancer','SSc-ILD','Systemic Sclerosis-Associated Interstitial Lung Disease','cHP','Chronic Hypersensitivity Pneumonitis','Pneumoconiosis','pneumoconiosis','Pulmonary Fibrosis','pulmonary fibrosis','pleural indentation','PI','Pulmonary nodules','Lung nodules')
#key2
data_t=c('scRNA-seq','scATAC-seq','snRNA-seq','singlecell RNA sequencing',
'spatial','stero-seq','visium',
'BTR','TCR'
)
all_query=list()
for(i in organs){all_query[[i]]=paste("(",i,")"," ","AND"," ","(",data_t,")",sep="")}
all_query1=unlist(all_query)


# loop for all 
combined_data <- data.frame()
for(i in 2:length(all_query1)){
  epm=epm_query(all_query1[i])
  epm <- epm_fetch(epm, format = 'xml')
  epm <- epm_parse(epm)
  proc_data <- get_epm_data(epm)
  combined_data <- rbind(combined_data, proc_data)
}
# save
write.csv(combined_data, "filename_combined_data.csv", row.names = FALSE)
