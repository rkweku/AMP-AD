# AMP-AD

##Description
This is an incomplete code which functions to identify a variety of patient  
samples within a Synapse project and organize them into one table with them  
into one table with the patient ID being a primary key. Only functions properly  
with the ROSMAP data for now. Takes clinical data and merges that with the key  
file, then it identifies the Synapse ID of various data files if it exists for  
that patient.

##Requirements
*Bio* - http://biopython.org/wiki/Main_Page  
*lxml.etree* - http://lxml.de/  
*pandas* - http://pandas.pydata.org/  
*synapseclient* - http://python-docs.synapse.org/

These may easily be installed using (Python) PIP. Intructions to install PIP -  
https://pip.pypa.io/en/stable/installing.html

##Notes
This code is not yet complete. It is necessary that the data for each study  
is parsed still. All files for a given project are currently identified and  
saved in a dictionary. The structure of the dictionary is a nested dictionary  
where the main dictionary has each study and each study has its own large  
dictionary with each unique data type as a key and a list of synapseIDs  
for all files of that data type. Currently, the script ignores all files that  
have 'None' as a data type

##Arguments
***synID:*** The Synapse ID of the project

##Contact
Reza Hammond  
RezaKweku@gmail.com
