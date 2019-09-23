# EQTL
Python Program (Library) to calculate linear regression between gene expression and gene synapse

# Dependencies
statsmodels, pyvcf, mysql-connector-python, scikit-learn, pandas, pysam

# Before Running: please install the following libraries
```{}
pip install --upgrade --no-deps statsmodels
```
```{}
pip install PyVCF
```
```{}
pip install mysql-connector-python
```
```{}
pip install scikit-learn
```
```{}
pip install pandas
```
```{}
pip install pysam
```

# Data

The Data file should be put into the Data folder

Acceptible data format:
*synpase*
.vcf, .txt, .vcf.gz

*gene expression*
.txt

**Important notes**
- Please put the tabix index file for vcf: *.vcf.gz.tbi*  inside the same directory.
- The name prefix of the tabix index file should be exactly same with the vcf file.

# Input
- User will nee to specify four input: 
    1. synapse name (id)
    2. gene expression name
    3. synapse vcf data path
    4. gene expression data path
- User will specify their input by modifying the *user_input.csv* inside the *user_input* folder.

# TODO
1. Optimize the indexing of gene expression file (currently use a map for indexing, higher space complexity).
2. Make sure the result is correcct when using real data to test
3. Handle more exceptions and corner cases. 