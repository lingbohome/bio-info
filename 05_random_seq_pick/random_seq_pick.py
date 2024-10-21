#! /usr/bin/env python3
#导入模块，初始传递命令
import argparse
import re
import os
import time
import random
#import pysam
parser = argparse.ArgumentParser(add_help = False, usage = '\npython3 bam_random_selection.py -i [inputbam] -s [sample_name] -n [read_num] -o [outpath] ')
required = parser.add_argument_group('必选项')
optional = parser.add_argument_group('可选项')
required.add_argument('-i', '--inputbam', metavar = '[input_path]', help = '输入文件路径', required = True)
required.add_argument('-s', '--sample_name', metavar = '[sample_name]', help = '样本名称', required = True)
required.add_argument('-n', '--read_num', metavar = '[read_num]', help = 'read数量', required = True)
required.add_argument('-o', '--outpath', metavar = '[output_path]', help = '输出文件路径', required = True)
optional.add_argument('-h', '--help', action = "help", help = "帮助信息")
args = parser.parse_args()
in_bamfile = args.inputbam 
n = args.read_num
s = args.sample_name
#def  RandomExtractReads(in_bamfile,n):
i = 0
samtools = os.getenv('SAMTOOLS_0_1_19')
out_bamfile_Header = args.outpath + "/" + s +".rmdup.target.pickseq.sam"
out_bamfile = args.outpath + "/" + s + ".rmdup.target.pickseq.txt"
out_newbam = args.outpath + "/" + s + ".rmdup.target.pickseq.bam"
f2 = open(out_bamfile_Header,'a+')
#print(out_bamfile_Header)
#print(out_bamfile)
os.system("%s view -H  %s > %s" %(samtools,in_bamfile,out_bamfile_Header))
os.system("%s view %s > %s" %(samtools,in_bamfile,out_bamfile))
#f2 = open(out_bamfile,'w+')
bam_list = []
with open(out_bamfile, 'r') as f:
	bam_list = f.readlines()
sam_num = len(bam_list)
f.close()
if (sam_num > int(n) ):
	sam_list = []
	sam_random_list = []
###产出list_id
	for i in range(sam_num):
		sam_list.append(i)
		#print(sam_list)
	sam_random_list = sorted(random.sample(sam_list,int(n)))
	#print(sam_random_list)
	for read_id in sam_random_list:
		#with open(in_bamfile) as bamfile:
		#	line_num =0
		#	for line in bamfile:
		#print(read_id)
		value = bam_list[read_id]
		value = value.strip()
		#print(value)
		#print(type(value))
		f2.write(value+'\n')
#print(bam_list)
	f2.close()
	del bam_list
	os.system("rm %s" %(out_bamfile))
	os.system("%s view -b -S %s > %s" %(samtools,out_bamfile_Header,out_newbam))
	os.system("rm %s" %(out_bamfile_Header))
	os.system("%s index %s" %(samtools,out_newbam))
else:
	os.system("cp %s %s " %(in_bamfile,out_newbam))
	os.system("rm %s %s" %(out_bamfile,out_bamfile_Header))
	os.system("%s index %s" %(samtools,out_newbam))
