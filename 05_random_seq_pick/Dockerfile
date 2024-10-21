FROM registry.cn-shanghai.aliyuncs.com/bio-cloud/bioinfo-aln:v0.2
WORKDIR /data1

ADD random_seq_pick.py random_seq_pick.sh /usr/local/bin/

# 赋予执行权限
RUN chmod +x /usr/local/bin/random_seq_pick.sh && chmod +x /usr/local/bin/random_seq_pick.py

# 设置工作目录
WORKDIR /data

