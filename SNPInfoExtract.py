## SNPInfoExtract.py 
## Extract SNP info from HMP file and pass the output to WhereTheSNP.py Script 
## Author: zhang fei 
## Please contact me: zhangfei-123@foxmail.com
## 2018-06-20

import sys

fl_snpquery = sys.argv[1]
fl_hmp = sys.argv[2]
#prefix = sys.argv[3]

## build snp_dic
with open(fl_hmp, 'r') as fh_hmp:
    dic_hmp = {}
    for line in fh_hmp:
        line = line.strip('\n').split('\t')
        line = line[0: 5]
        if dic_hmp.get(line[0]) is None:
            dic_hmp[line[0]] = [line[0], line[2], line[3], line[1], line[4]]

with open(fl_snpquery, 'r') as fh_snpquery:
    result_snpquery = []
    for line in fh_snpquery:
        line = line.strip('\r\n')
        if dic_hmp.get(line) is not None:
            result_snpquery.append(dic_hmp[line])
            # print line, " is in the hmp file."
        else:
            print line, " not in the hmp file, please check the hmp or input file."

with open("./4WhereTheSNP.txt", 'w') as fh_output:
    result_snpquery = ['\t'.join(i) for i in result_snpquery]
    fh_output.write('\r\n'.join(result_snpquery))
