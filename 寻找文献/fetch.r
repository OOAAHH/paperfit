library(easyPubMed)
#Use R env R432
#key1
organs=c('Lung','nose','pharynx','larynx','trachea','bronchi','alveolar','COVID-19','SARS','flu',
        
        'Blood',
        'Nose',
        'Mouth',
        'Pharynx',
        'Larynx',
        'Trachea',
        'Bronchi',
        'Lungs',
        'Lung',
        'Alveoli',#新增
        'Diaphragm',#新增
        'Respiratory system',#新增'respiratory system',#新增

        'Asthma',
        'Bronchitis',
        'Chronic Obstructive Pulmonary Disease',
        'Pneumonia',
        'Tuberculosis',
        'Lung Cancer',
        'Pulmonary Fibrosis',
        'Pulmonary Hypertension',
        'Cystic Fibrosis',
        'Acute Respiratory Distress Syndrome',
        'Pleural Effusion',
        'Pneumothorax',
        'Sarcoidosis',
        'Lung Abscess',
		'Flu',
		'Influenza',
		'COVID-19',
		'SARS-cov-2',
		'SARS',
		'Allergy',
		'Asthma',
		'Silicosis',
		'Pneumoconiosis',
		'Pulmonary Fibrosis',

        'COPD','Chronic Obstructive Pulmonary Disease',
        'ILD','Interstitial Lung Disease',
        'IPF','Idiopathic Pulmonary Fibrosis',
        'KD','Kawasaki Disease',
        'LSCC','Lung Squamous Cell Carcinoma',
        'LUAD','Lung Adenocarcinoma',
        'NSCLC','Non-Small Cell Lung Cancer',
        'NSIP','Nonspecific Interstitial Pneumonia',
        'SCC','Squamous Cell Carcinomas',
        'SCLC','Small Cell Lung Cancer',
        'SSc-ILD','Systemic Sclerosis-Associated Interstitial Lung Disease',
        'cHP','Chronic Hypersensitivity Pneumonitis',
        )
#key2
data t=c('scRNA-seq','scATAC-seq','snRNA-seq','singlecell RNA sequencing'
        'spatial','stero-seq','visium',
        'BTR','TCR',
        )
# 保存organs_PM_df为CSV文件
write.csv(organs_PM_df, "/home/sunhao/RAW/filename.csv", row.names = FALSE)
