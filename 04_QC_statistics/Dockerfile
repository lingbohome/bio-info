FROM registry.cn-shanghai.aliyuncs.com/bio-cloud/bioinfo-aln:v0.2
WORKDIR /data1

RUN wget https://github.com/arq5x/bedtools2/releases/download/v2.25.0/bedtools-2.25.0.tar.gz && \
    tar -zxvf bedtools-2.25.0.tar.gz && \
    rm bedtools-2.25.0.tar.gz && \
    cd bedtools2 && \
    cp /usr/bin/python3 /usr/bin/python && \
    make && \
    mkdir -p $SOFT/bedtools2 && cp -r bin/ $SOFT/bedtools2/ && \
    cd ../ && rm -rf bedtools2/

RUN cpan Getopt::Long::Descriptive    

ENV Bedtools $SOFT/bedtools2/bin/bedtools 

ADD library_qc_ctDNA.pl QC_statistics.sh /usr/local/bin/

# 赋予执行权限
RUN chmod +x /usr/local/bin/QC_statistics.sh && chmod +x /usr/local/bin/library_qc_ctDNA.pl

# 设置工作目录
WORKDIR /data

