"""
This script has been written to export legacy Array data from moka 
"""
import os 
import pandas as pd
from ConfigParser import ConfigParser 
import pyodbc 

# Read config file(must be called config.ini and stored in the same directory as script)
config_parser = ConfigParser()
print_config = config_parser.read(os.path.join(os.path.dirname(os.path.realpath(__file__)), "config.ini"))

# Create pyodbc connection to moka 
cnxn = pyodbc.connect('DRIVER={{SQL Server}}; SERVER={server}; DATABASE={database};'.format(
         server=config_parser.get("MOKA", "SERVER"),
         database=config_parser.get("MOKA", "DATABASE")))
         

'''=================================================== SCRIPT =================================================== ''' 
# Query to pull data from moka 
export_patient_data_SQL = ("SELECT DISTINCT 'chr' + [Chr] AS Chromo, ArrayOligoPreliminaryResults.Start19," 
    "ArrayOligoPreliminaryResults.Stop19, Patients.PatientID, ArrayOligoPreliminaryResults.Pathogenic, "
    "Referral.Referral, Status.Status, ArrayOligoPreliminaryResults.CNVTypeID, "
    "[dbo].[gwv-patientlinked].[Gender], Patients.Sexed, Patients.BookinSex, " # Need all three to get all patient's sex
   " 'Patient result for' + ' ' + [Patients].[PatientID] AS PatientResult, "
    "[band19] + '(' + CONVERT(varchar(20), start19) + '-' + CONVERT(varchar(20), stop19) + ')' + [Change] AS [DESCtwo], " # Convert integers to make a string
    "Arrays.DateRecieved, [Patients].[PatientID] + [AltForBED] AS PatientIDInheritance, "
    "Phenotype, ArrayOligoPreliminaryResults.Copies, Change.Change, ArrayOligoPreliminaryResults.WholeChromosome, "
    "Patients.BookinDOB,  [dbo].[gwv-patientlinked].[DoB], ArrayTest.RequestedDate  " # To calculate patient's age 
    "FROM (((((((((((Patients INNER JOIN ResultCode ON Patients.OverallResultCodeID = ResultCode.ResultCodeID) "
    "INNER JOIN ArrayOligoPreliminaryResults ON Patients.InternalPatientID = ArrayOligoPreliminaryResults.InternalPatientID) "
    "INNER JOIN Chromosome ON ArrayOligoPreliminaryResults.ChrID19 = Chromosome.ChrID) " 
    "INNER JOIN ArrayLabelling ON ArrayOligoPreliminaryResults.DNALabellingID = ArrayLabelling.DNALabellingID) "
    "INNER JOIN Arrays ON ArrayLabelling.ArrayID = Arrays.ArrayID) " 
    "INNER JOIN Status ON ArrayOligoPreliminaryResults.Pathogenic = Status.StatusID) "
    "INNER JOIN [dbo].[gwv-patientlinked] ON [dbo].[gwv-patientlinked].[PatientTrustID] = Patients.PatientID) " 
    "INNER JOIN Change ON Change.ChangeID = ArrayOligoPreliminaryResults.Copies) " # one patient has no copy information
    "INNER JOIN Phenotype ON Phenotype.InternalPatientID = ArrayOligoPreliminaryResults.InternalPatientID) " # Only return patients with phenotype information
    "INNER JOIN ArrayTest ON ArrayTest.InternalPatientID = Patients.InternalPatientID)"
    "INNER JOIN Referral ON ArrayTest.ReferralID = Referral.ReferralID)"
    "LEFT JOIN Inheritance ON Inheritance.InheritanceID = ArrayOligoPreliminaryResults.InheritanceID " # Some patients have no inheritance reported
    "WHERE (((Referral.ReferralID)=1199901171 OR (Referral.ReferralID)=2 OR (Referral.ReferralID)=1199901185 " # Only referrals for 'POC/Tissues', 'General Ref', 'Prenatal' Or
    "OR (Referral.ReferralID)=1199901199 OR (Referral.ReferralID)=1199901202 OR (Referral.ReferralID)=1199901195 " # 'Repeat POC', 'repeat prenatal', 'repeat gen ref' Or
    "OR (Referral.ReferralID)=1199901215) " # referred as a 'pseudorush'
    "AND ((ArrayOligoPreliminaryResults.CNVTypeID)=1190384964 OR (ArrayOligoPreliminaryResults.CNVTypeID)=1190384990 " # Only 'reported', 'Not reported (unrelated to referral)' Or
    "OR (ArrayOligoPreliminaryResults.CNVTypeID)=1190384994) " #  'Not reported (insufficient evidence)'
    "AND ((Patients.OverallResultCodeID)<>1189679594 And (Patients.OverallResultCodeID)<>1189679593)  " # Results code is not 'Preliminary result' or 'Fail / testing not possible'
    "AND ((ArrayOligoPreliminaryResults.Pathogenic)=1202218781 OR " # Patient result is 'Uncertain clinical significance (class 3)' Or
    "(ArrayOligoPreliminaryResults.Pathogenic)=1202218783 Or (ArrayOligoPreliminaryResults.Pathogenic)=1202218788) " # 'Likely to be pathogenic (class 4)' Or 'Pathogenic / abnormal result (class 5)'
    "AND ((Arrays.DateRecieved) BETWEEN '2016-01-01 00:00:00' AND '2021-12-31 00:00:00') " # For current array plaform
    "AND ((ArrayTest.StatusID) = 4)" # Only tests which are completed 
    "AND ((Patients.BookinDOB)IS NOT NULL OR ([dbo].[gwv-patientlinked].[DoB]) IS NOT NULL))") # Only patients with a DOB in either Moka or GW
# Run SQL query and store results in a pandas dataframe 
export_patient_data_df = pd.read_sql(export_patient_data_SQL, cnxn)
print("Data pulled from moka")

# Save as a tsv file, some characters are not ascii so encoding in utf-8
export_patient_data_df.to_csv("moka_array_export_220215.tsv",sep = "\t", encoding='utf-8')
print("Data saved")
