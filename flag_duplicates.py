
# coding: utf-8

# In[22]:

#Imports & Setup

"""
Description: Flags duplicate patients in the argument .csv file by adding a column called 'Flag' that maps to
             unique user tests. Thus, the number in the 'Flag' column for any given row will only match to duplicate
             persons

Sample Usage: Use 'python3 flag_duplicates.py [$filename$.csv]' to flag all duplicates and output a new .csv to the path
./filename.flagged_duplicates.csv, where filename is replaced with a 
"""

import pandas as pd
import numpy as np
import math as mt
import sys
from fuzzywuzzy import fuzz
from fuzzywuzzy import process


# In[23]:

#Loading in data

#Ensuring that a file path has been passed into the script
improper_arg_msg = "Improper arguments passed in - use 'python3 flag_duplicates.py [$filename$.csv]'"
if len(sys.argv) < 2:
    assert False, improper_arg_msg
elif ".csv" not in sys.argv[1]:
    print("Using the sample dataset...\n", improper_arg_msg)
    sys.argv[1] = "deduplicator_sample_data_scramble.csv"


#Getting the roster and homework response paths
path = sys.argv[1]

lab_confirmed_flu = pd.read_csv(sys.argv[1])
lab_confirmed_flu


# In[24]:

#Find Fuzzy string matches for each Patient
try:
    patients = lab_confirmed_flu["Patient"]
except:
    lab_confirmed_flu["Patient"] = np.core.defchararray.add(lab_confirmed_flu['last_name'], lab_confirmed_flu['first_name'])
    patients = lab_confirmed_flu["Patient"]

#The first match of each list will be the row the patient was from
matches = [process.extract( query=patient, choices=patients, limit=max(25, int(len(patients) ** .5)) )
           for patient in patients]
matches


# In[28]:

#Find true duplicates and match them up

#Maps frozenset tuples of 2 filtered match indices >> patient index (-1 if they are not true matches)
filtered_match_ids = dict()
#The column to be added, containing the updated patient_id for each index
flags = np.arange( len(patients) )
#Ensures that one patient's index is not set multiple times (a single patient's index should not match via manual and automated detection more than once)
already_matched = set()
#Keeps track of all matches for displaying later
all_matches_in_dataset = []

def validate_matches(filtered_match_index, filtered_match_ids):
    """For a given filtered_match_index, returns whether all keys in filtered_match_ids that have filtered_match_index within the key have the same value
        :param int filtered_match_index : the index of a match
        :param dict filtered_match_ids  : Maps frozenset tuples of 2 filtered match indices >> patient index (-1 if they are not true matches)
        :return boolean                 : whether all keys in filtered_match_ids that have filtered_match_index within the key map to the same value
    """
    match_keys = [key for key in filtered_match_ids if filtered_match_index in key and filtered_match_ids[key] != -1]
    reference_value = filtered_match_ids[ match_keys[0] ]
    for key in match_keys:
        if filtered_match_ids[key] != reference_value: return False  
    return True
    
for patient_match_list in matches:
    patient = patient_match_list[0]
    patient_index = int(patient[2])
    
    #Narrow down all fuzzy string scores to only potential duplicates of "patient"
    all_matches = np.asarray(patient_match_list[1:])
    filtered_matches = all_matches[ np.asarray([int(match[1]) > 65 for match in all_matches]) ]
    
    #For all filtered matches, find true duplicates and give them the same patient_id
    #Note: simply because another patient passed the filter does NOT mean they are a true match 
    for filtered_match in filtered_matches:
        
        filtered_match_index = int(filtered_match[2])
        possible_match_key = frozenset([patient_index, filtered_match_index])
        
        #If this possible_match_key has already been checked and IS NOT a match
        if possible_match_key in filtered_match_ids and filtered_match_ids[possible_match_key] == -1:

            continue
            
        #If this possible_match_key has already been checked and IS a match, use the same value
        elif possible_match_key in filtered_match_ids:
            
            flags[patient_index] = filtered_match_ids[possible_match_key]
            break
        
        #If the possible_match_key has not already been checked, determine whether it is a true match
        else:
            patient_row        = lab_confirmed_flu.iloc[patient_index]
            filtered_match_row = lab_confirmed_flu.iloc[filtered_match_index]
            
            #Uncertain based on data -- ask user to take a closer look
            if ( (type(patient_row.DOB)               == float and np.isnan(patient_row.DOB ))        or
                 (type(filtered_match_row.DOB)        == float and np.isnan(filtered_match_row.DOB )) or
                 (type(patient_row.Collected)         == float and np.isnan(patient_row.Collected  )) or   
                 (type(filtered_match_row.Collected)  == float and np.isnan(filtered_match_row.Collected))
               ):
                
                msg = """Please press 'Y' if the two patients are matches and anything else if they are not: """
                print("\n\n===========================================================\nPlease examine the following:\n")
                print("\tTarget Patient:\n", patient_row)
                print("\tPotential match:\n", filtered_match_row)
                is_match = input("\n"+msg).strip().lower() == 'y'
                
            #Highly probable matches based on DOB + Collection Time + Test
            else:                

                DOB_match       = patient_row["DOB"]       == filtered_match_row["DOB"]
                Collected_match = patient_row["Collected"] == filtered_match_row["Collected"]
                
                try:
                    Test_match  = patient_row["Test"]      == filtered_match_row["Test"]
                
                except: #cc.dedup
                    try:
                        Test_match  = patient_row['Result']== filtered_match_row['Result']   
                    
                    except: #cho.a.dedup
                        try:
                            Test_match  = patient_row['flua']  == filtered_match_row['flua'] and patient_row['flub']  == filtered_match_row['flub'] 
                        
                        except: #cho.b.dedup
                            try:
                                Test_match  = (patient_row['influenza.a.h1']  == filtered_match_row['influenza.a.h1'] and 
                                               patient_row['influenza.a.h3']  == filtered_match_row['influenza.a.h3'] and 
                                               patient_row['x2009.inf.a.h1n1.rvp']  == filtered_match_row['x2009.inf.a.h1n1.rvp'] and
                                               patient_row['flu.b']  == filtered_match_row['flu.b'] and
                                               patient_row['rsv.a']  == filtered_match_row['rsv.a'] and
                                               patient_row['rsv.b']  == filtered_match_row['rsv.b'] and
                                               patient_row['parainfluenza.1']  == filtered_match_row['parainfluenza.1'] and
                                               patient_row['parainfluenza.2']  == filtered_match_row['parainfluenza.2'] and
                                               patient_row['parainfluenza.3']  == filtered_match_row['parainfluenza.3'] and
                                               patient_row['rhinovirus']  == filtered_match_row['rhinovirus'] and
                                               patient_row['adenovirus']  == filtered_match_row['adenovirus'] and
                                               patient_row['metapneumovirus']  == filtered_match_row['metapneumovirus']
                                              )
                            except:
                                assert False, "All cases should've been covered"
                                print(path)
                                print("\n\nSample Patient Row:\n\n")
                                print(patient_row)
                            
                is_match = DOB_match and Collected_match and Test_match

            #Storing the result of our comparison in filtered_match_ids
            if is_match:
                
                all_matches_in_dataset.append(patient_row)
                
                contradiction_msg = "The newly matched patient has already been matched -- this is a contradiction. Patient: " + str(patient)
                assert (filtered_match_index not in already_matched) or validate_matches(filtered_match_index, filtered_match_ids), contradiction_msg 
                
                filtered_match_ids[possible_match_key] = flags[patient_index]
                flags[filtered_match_index] = flags[patient_index]
                already_matched.add(filtered_match_index)
            
            else:           
                filtered_match_ids[possible_match_key] = -1
        
lab_confirmed_flu["Flag"] = flags 
lab_confirmed_flu


# In[29]:

#Saving the result

assert path[-4:] == ".csv"
new_path = path[:-4] + ".flagged_duplicates.csv"

lab_confirmed_flu.to_csv(path_or_buf=new_path, index=False)


# In[30]:

if len(all_matches_in_dataset) > 0:
    print("Dataset:", new_path)
    print("All matches in dataset:\n\n", all_matches_in_dataset)


# In[ ]:



