"""
*******************************************************************************************************
VERSION V1.0
1. Update tbi reading
2. !!!!Only accept file input => users would input their snp_id, snp_file_path, gene_exp_name, geneexp_file_path
    inside the user_input folder => modifying user_input.csv to run the program
3. 
*******************************************************************************************************
"""

# INSTRUCTIONS
# install dependencies: statsmodels, pyvcf, mysql-connector-python, scikit-learn, pandas, pysam
# 
# pip install --upgrade --no-deps statsmodels
# pip install PyVCF
# pip install mysql-connector-python
# pip install scikit-learn
# pip install pandas
# pip install pysam
# 

# TODO: Optimize the indexing of gene_exp table, currently build a map for indexing => might take large space

from sklearn import linear_model
import pandas as pd
import statsmodels.formula.api as sm
import mysql.connector
import vcf
import pysam
import csv



# --------------------------------------------------------------------------------------------------------------
# FUNCTIONS for parsing vcf format data file:

# get corresponding gt data for input snp_id
# input a snp_id and file_path
# return a dictionary of sample_name: gt_value
def process_data(snp_id, file_path):
    # Step 1: map snp_id to chrom name and snp_position
    query = id_to_pos(snp_id)
    start_pos = int(query[1])
    end_pos = int(query[2])
    chrom = query[0]
    # print(chrom,start_pos,end_pos)
    # Step 2: read vcf file
    vcf_reader = vcf.Reader(filename=file_path)
    # print(vcf_reader)    

    # Step 3: search for corresponding snp_id in the vcf file using tbi and return it record
    try:
        record = next(vcf_reader.fetch(chrom=chrom[3:],start=start_pos-1,end=end_pos-1))
        result = process_vcf(record)
        return result
    except:
        print("Search for input snp_id failed!!")
        print("Please make sure tbi file and vcf file are in the same folder and check the path or name tbi file position. If no error in path and names present, try remake the tbi file.")
        print("------------------------------------------------------------")
        print("Do exaustive search instead...... press Ctrl C to force stop.")

        for record in vcf_reader:
            if record.ID == snp_id:
                result = process_vcf(record)
                return result
        
        print("snp_id cannot be found in the file, check the input file for correctness.")



# helper function to map snp_id to chrom name and snp_postion
def id_to_pos(snp_id):
    return query_db(connect_db(), snp_id)


# helper function to find samples' gt_value from chrom name and snp_position:
def process_vcf(record):
    result = {}
    ref = record.REF
    # print(ref)
    alt = record.ALT
    # print(alt)
    for sample in record.samples:
        gt = sample['GT']
        if gt[0] != '.' and gt[2] != '.':
            if int(gt[0]) < 2 and int(gt[2]) < 2:
                value_sum = int(gt[0]) + int(gt[2])
                result[sample.sample] = value_sum
    return [result, [ref, alt]]


# helper function to connect to gene browser database
def connect_db():
    try:
        db = mysql.connector.connect(user="genome",
                                     host="genome-mysql.soe.ucsc.edu",
                                     database="hg38",
                                     port="3306")
        return db
    except:
        print("connection failed, try another port")


# Helper function to query snp_position with snp_id from the database
def query_db(db, snp_id):
    # use cursor() method to find the coordinate

    cursor = db.cursor()

    # SQL query
    snp_id = '("' + snp_id + '")'
    sql = 'SELECT chrom, chromStart, chromEnd, name FROM snp151 where name in ' + snp_id
    res = []
    try:
        # execute sql statement
        cursor.execute(sql)
        # fetch chrom, start and end positions
        results = cursor.fetchall()
        for result in results:
            chrom = result[0]
            start_pos = result[1]
            end_pos = result[2]
            rs_id = result[3]
            res.append(chrom)
            res.append(start_pos)
            res.append(end_pos)
            res.append(rs_id)
        cursor.close()
        db.close()
        return res
    except:
        print("Error: unable to fecth data")



# --------------------------------------------------------------------------------------------------------------
# FUNCTIONS for parsing gene expressions data file

# Build HashMap for expression input file:
def format_data(file_name):
    try:
        data = open(file_name, 'r')
        rows = data.readlines()
        result = {}
        count = 0
        for row in rows:
            if row[:4] != "Name":
                count += 1
            else:
                break
        header = rows[count]
        header_list = header.split()
        sample_names = header_list[2:]
        body = rows[count + 1:]
        for row in body:
            row_list = row.split()
            if row:
                result[row_list[1]] = row_list[2:]
        data.close()
        return [result, sample_names]
    except FileNotFoundError:
        print("File not found! Please input the correct file path")
        return


# This is the function to find the sample data using gene_expression name
def find_gene_exp(gene_exp, file_name):
    result = {}
    # Load data in
    try:
        data = format_data(file_name)[0]
        sample_names = format_data(file_name)[1]
    except FileNotFoundError:
        return

    # build the result dictionary for pandas data frame
    #  key: sample name -> value: sample data
    try:
        # query the gene expression name from the index hash map (dictionary)
        exp_data = data[gene_exp]
        for i in range(len(sample_names)):
            result[change_name(sample_names[i])] = int(exp_data[i])
    except KeyError:
        print("Gene expression not found! Please Try Again")
        return
    return result


# helper function to change the sample to name to match the ones in vcf files
def change_name(name):
    result = name.split('-')
    return result[0] + "-" + result[1]



# --------------------------------------------------------------------------------------------------------------
# FUNCTIONS for merging the two data into a same file


# merging data
# params: snp dict, gene_exp dict: key => name, value => values
def merge_data(snp, gene_exp):
    # df_snp = pd.DataFrame([(item[0],item[1][0]) for item in snp.items()], columns=['sample_name', 'snp_value'])
    df_snp = pd.DataFrame(list(snp.items()), columns=['sample_name', 'snp_value'])
    gene_exp_df = pd.DataFrame(list(gene_exp.items()), columns=['sample_name', 'gene_exp_value'])
    result = pd.merge(df_snp, gene_exp_df, how='inner', on='sample_name')
    return result



# --------------------------------------------------------------------------------------------------------------
# MAIN PROGRAM starts here

def main():

    # Ask the user for inputs
    # 1. snp_name 2. snp_vcf file path 3.
    # snp_name = input("please input the snp name (snp_RSid)")
    # print(snp_name)
    # gene_exp_name = input("please input gene expression name")
    # print(gene_exp_name)
    # snp_file_path = input("please input snp data path (Accept compressed or not compressed vcf format; "
    #                       "for faster search, please put the tbi index file in the same directory as the vcf file)")
    # print(snp_file_path)
    # gene_exp_file_path = input("please input gene expression data path (Accept txt file input)")
    # print(gene_exp_file_path)
    # parameter_file = input("please input the parameter file")
    input_file = open("user_input/user_input.csv")
    parameter_list=input_file.readlines()[1].split(',')
    print("------------------------------------------------------------")
    print("Show loaded parameter list:")
    print(parameter_list)
    # for test
    snp_name = parameter_list[0]
    gene_exp_name = parameter_list[1]
    snp_file_path = parameter_list[2]
    gene_exp_file_path = parameter_list[3]

    # Loading Data
    print("Loading Data...")
    print("------------------------------------------------------------")
    snp_dict = process_data(snp_name, snp_file_path)
    gene_exp_dict = find_gene_exp(gene_exp_name, gene_exp_file_path)
    snp_ref = snp_dict[1][0]
    snp_alt = snp_dict[1][1]
    snp_dict = snp_dict[0]

    # Merge two dicts based on sample name
    df = merge_data(snp_dict, gene_exp_dict)

    # Building model
    print("Building Model...")
    print("------------------------------------------------------------")

    # Using sklearn
    # model = linear_model.LinearRegression()
    # X = df[['gene_exp_value']]
    # y = df['snp_value']
    # model.fit(X, y)
    # print("coefficients: " + str(model.coef_))
    # print("intercept: " + str(model.intercept_))
    # print("R^2" + str(model.score(X, y)))

    # stats_model
    res = sm.ols('snp_value ~ gene_exp_value', data=df).fit()
    # print(res.params)
    # print(res.pvalues)

    # Build result
    print("Result:")
    result = {"snp_id": snp_name, "snp_REF": snp_ref, "snp_ALT": snp_alt, "gene_expression": gene_exp_name,
              "intercept": res.params[0], "intercept_pvalue": res.pvalues[0], "beta": res.params[1],
              "beta_pvalue": res.pvalues[1]}
    print(result)


if __name__ == '__main__':
    main()

# --------------------------------------------------------------------------------------------------------------
# TEST CASES

# --------------------------------------------------------------------------------------------------------------
# 1. Tests for parsing gene expression data file

# print(find_gene_exp("MIR1302-", "Data/gene_exp.txt"))
# print(format_data("Data/gene_exp.txt"))
# data = format_data("Data/gene_exp.txt")
# exp_data = data[0]
# print(exp_data["a"])

# --------------------------------------------------------------------------------------------------------------
# 2. Tests for parsing vcf format data file

# 1182961
# print(process_data("rs6040355", "Data/example.vcf.gz"))
# print(process_data("rs607852", "Data/ALL.chr1.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz"))
# process_data("rs6040355", "Data/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
# test = vcf.Reader(filename="Data/ALL.chr20.phase3_shapeit2_mvncall_integrated_v5a.20130502.genotypes.vcf.gz")
# record = test.fetch('20', 1182960, 1182961) # doctest: +SKIP
# record=next(record)
# print(process_data(record)))

# --------------------------------------------------------------------------------------------------------------
# 4. Tests for mergind two sets of data

# snp = vcf_data.process_data("rs6040355", 'Data/example.vcf')
# gene_exp = gtex_input.find_gene_exp("MIR1302-11", "Data/gene_exp.txt")
# print(snp)
# print(gene_exp)
#
# df = merge_data(snp, gene_exp)
# print(df)
# ex_record = vcf_input.find_record("rs6054257", "Data/vcf_data.txt")
# snp = vcf_input.manipulate_record(ex_record)
# gene_exp = gtex_input.find_gene_exp("MIR1302-11", "Data/gene_exp.txt")

# df = merge_data(snp, gene_exp)
# print(df)
# print(df["snp_value"])
# print(df["gene_exp_value"])
