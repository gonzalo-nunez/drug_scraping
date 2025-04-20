from fastapi import FastAPI, Query
from fastapi.middleware.cors import CORSMiddleware
from utils import *
import pandas as pd

app = FastAPI()

app.add_middleware(
    CORSMiddleware,
    allow_origins=['*'],
    allow_credentials=True,
    allow_methods=['*'],
    allow_headers=['*'],
)

@app.get('/clinical-data')
def get_clinical_data(condition: str):
    df, total = gather_data(condition)
    df = df.fillna('')
    status, interventions, drug_list = gather_stats(df)
    #drug_df = obtain_drug_data(extract_drugs(df))
    #drug_df = drug_df.fillna('')
    
    return {
        'total_studies': total,
        'studies_df': df.to_dict(orient='records'),
        'status_stats': status,
        'interventions_stats': interventions,
        'drug_stats': drug_list#,
        #'drug_data': drug_df.to_dict(orient='records')
    }