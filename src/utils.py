import requests
import pandas as pd
import matplotlib.pyplot as plt
import numpy as np
import re
import difflib
import random
from bs4 import BeautifulSoup
import requests
import pandas as pd
import ast
from io import StringIO
from collections import Counter
from chembl_webresource_client.new_client import new_client

def clean_text(text):
  new_text = re.sub(r"[\['\]]", "", text)
  return new_text

def clean_df(df):
  for col in df.columns:
    df[col] = df[col].astype(str).apply(clean_text)
  df = df[df['Interventions'].str.contains('DRUG') & (df['Study Type']=='INTERVENTIONAL') & (df['Phases']=='PHASE3')]
  return df

def gather_data(condition):
    req = requests.get(f'https://clinicaltrials.gov/api/v2/studies?format=json&query.cond={condition}&fields=NCTId&countTotal=true&pageSize=1')
    data = req.json()
    StudiesFound = data['totalCount']
    print(StudiesFound,f'studies found using the "{condition}" search parameter')
    req = requests.get(f'https://clinicaltrials.gov/api/v2/studies?format=csv&query.cond={condition}&fields=NCT+Number%7CConditions%7CStudy+Title%7CBrief+Summary%7CStudy+Type%7CInterventions%7CPhases%7CStudy+Status%7CCompletion+Date%7CSponsor&pageSize=1000')
    csv_data = StringIO(req.text)
    df = pd.read_csv(csv_data)
    while(True):
        try:
            next_page_token = req.headers['x-next-page-token']
            req = requests.get(f'https://clinicaltrials.gov/api/v2/studies?format=csv&query.cond={condition}&fields=NCT+Number%7CConditions%7CStudy+Title%7CBrief+Summary%7CStudy+Type%7CInterventions%7CPhases%7CStudy+Status%7CCompletion+Date%7CSponsor&pageSize=1000&pageToken={next_page_token}')
            csv_data = StringIO(req.text)
            df2 = pd.read_csv(csv_data, header= None)
            df2.columns = df.columns
            df = pd.concat([df, df2], axis=0, ignore_index=True)
        except:
            df = df.reset_index(drop=True)
            break
    return df, StudiesFound

def gather_stats(df):
    total_values = len(df)
    temp_df = pd.DataFrame()
    #Study Status (Completed, recruiting, etc)
    status_count = df['Study Status'].value_counts()
    remaining_count = total_values - int(sum(status_count.values))
    status = [status_count.index.tolist()+['Not Reported'], status_count.values.tolist()+[remaining_count]]

    #Intervention types (Drug, Device, etc)
    temp_df['Interventions'] = df['Interventions'].str.split(':').str[0]
    interv_count = temp_df['Interventions'].value_counts()
    remaining_count = total_values - int(sum(interv_count.values))
    interventions = [interv_count.index.tolist()+['Not Reported'], interv_count.values.tolist()+[remaining_count]]

    #Drug names
    drugs = []
    for field in df[df['Interventions'].str.contains('DRUG', na=False)]['Interventions'].tolist():
        value = field.split('|')
        for a in value:
            b = a.split(': ')
            if b[0] == 'DRUG':
                drugs.append(b[1])
    drugs = [drug for drug in drugs if 'placebo' not in drug.lower()]
    counts = Counter(drugs).most_common()
    labels, values = zip(*counts)
    drug_list = [list(labels), list(values)]
    return status, interventions, drug_list

def extract_drugs(df):
    NCTId = []
    Condition = []
    DrugName = []
    StudyStatus = []
    CompletionYear = []
    for i in df.index:
        try:
            Interventions = df['Interventions'][i].split('|')
            for intervention in Interventions:
                a = intervention.split(': ')
                if a[0] == 'DRUG' and 'placebo' not in a[1].lower() and 'saline' not in a[1].lower():
                    NCTId.append(df['NCT Number'][i])
                    Condition.append(df['Conditions'][i])
                    DrugName.append(a[1])
                    StudyStatus.append(df['Study Status'][i])
                    year = df['Completion Date'][i]
                    CompletionYear.append(year.split('-')[0] if pd.notna(year) else np.nan)
        except:
            _
    new_data = {'Drug':DrugName, 'Condition':Condition, 'NCTId':NCTId, 'Study Status':StudyStatus, 'Completion Year':CompletionYear}
    new_df = pd.DataFrame(new_data)
    return new_df

def obtain_drug_data(df):
    NCT = []
    drug = []
    approval = []
    ChID = []
    struct = []
    moltype = []
    names = []
    df = df.drop_duplicates(subset='Drug', keep='first')
    df = df.reset_index(drop=True)
    molecule = new_client.molecule
    for i in df.index:
        Drugs = df['Drug'][i].split(' ')
        for d in Drugs:
            mols = molecule.filter(molecule_synonyms__molecule_synonym__iexact=d).only(['molecule_chembl_id', 'molecule_structures', 'molecule_type', 'first_approval__gte', 'pref_name'])
            if len(mols)>0:
              m = mols[0]
              drug.append(d)
              approval.append(int(m['first_approval']) if pd.notna(m['first_approval']) else np.nan)
              ChID.append(m['molecule_chembl_id'])
              struct.append(m['molecule_structures']['canonical_smiles'] if pd.notna(m['molecule_structures']) else np.nan)
              moltype.append(m['molecule_type'])
              names.append(m['pref_name'])
              NCT.append(df['NCTId'][i])
    dlist = {'Drug':drug, 'Approval Date':approval, 'ChEMBL ID':ChID, 'Molecule Type':moltype, 'Structure':struct, 'Prefered Name':names, 'NCTId':NCT}
    drug_list = pd.DataFrame(dlist)
    return drug_list

