FROM registry.cn-shanghai.aliyuncs.com/bio-cloud/bioinfo-aln:v0.2
WORKDIR /data1

ENV MSISENSOR $SOFT/msisensor/msisensor
ENV MSISENSOR_PRO $SOFT/msisensor-pro/msisensor-pro

RUN mkdir -p $SOFT/msisensor && wget https://github.com/ding-lab/msisensor/releases/download/0.6/msisensor.linux -O $MSISENSOR && chmod a+x $MSISENSOR
RUN apt-get update && apt-get install -y libhts-dev
RUN mkdir -p $SOFT/msisensor-pro && wget https://github.com/xjtu-omics/msisensor-pro/releases/download/v1.3.0/msisensor-pro-v1.3.0 -O $MSISENSOR_PRO && chmod a+x $MSISENSOR_PRO
##ENV LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu/libhts.so    

##RUN apt-get update && apt-get install -y libcurl4-openssl-dev libncursesw5-dev
# RUN wget https://github.com/xjtu-omics/msisensor-pro/archive/refs/tags/v1.3.0.tar.gz -O msisensor-pro-1.3.0.tar.gz && \
#     tar -xzvf msisensor-pro-1.3.0.tar.gz && \
#     rm msisensor-pro-1.3.0.tar.gz && \
#     cd msisensor-pro-1.3.0 && \
#     ./INSTALL && mkdir -p $SOFT/msisensor-pro && mv binary/msisensor-pro $MSISENSOR_PRO && \
#     cd ../ && rm -rf msisensor-pro-1.3.0/  

# 设置工作目录
WORKDIR /data

