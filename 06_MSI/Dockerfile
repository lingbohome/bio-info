FROM registry.cn-shanghai.aliyuncs.com/bio-cloud/bioinfo-aln:v0.2
WORKDIR /data1

RUN apt-get update && apt-get install -y python3-pip libhts-dev bc libstdc++6 &&  rm -rf /var/lib/apt/lists/*

ENV MSISENSOR $SOFT/msisensor/msisensor
ENV MSISENSOR_PRO $SOFT/msisensor-pro/msisensor-pro

RUN mkdir -p $SOFT/msisensor && wget https://github.com/ding-lab/msisensor/releases/download/0.6/msisensor.linux -O $MSISENSOR && chmod a+x $MSISENSOR

RUN mkdir -p $SOFT/msisensor-pro && wget https://github.com/xjtu-omics/msisensor-pro/releases/download/v1.3.0/msisensor-pro-v1.3.0 -O $MSISENSOR_PRO && chmod a+x $MSISENSOR_PRO
##ENV LD_LIBRARY_PATH /usr/lib/x86_64-linux-gnu/libhts.so::$LD_LIBRARY_PATH    

##RUN apt-get update && apt-get install -y libcurl4-openssl-dev libncursesw5-dev
# RUN wget https://github.com/xjtu-omics/msisensor-pro/archive/refs/tags/v1.3.0.tar.gz -O msisensor-pro-1.3.0.tar.gz && \
#     tar -xzvf msisensor-pro-1.3.0.tar.gz && \
#     rm msisensor-pro-1.3.0.tar.gz && \
#     cd msisensor-pro-1.3.0 && \
#     ./INSTALL && mkdir -p $SOFT/msisensor-pro && mv binary/msisensor-pro $MSISENSOR_PRO && \
#     cd ../ && rm -rf msisensor-pro-1.3.0/



# ENV CONDA_DIR $SOFT/anaconda3
# RUN wget --quiet https://repo.anaconda.com/archive/Anaconda3-2023.09-0-Linux-x86_64.sh -O ~/anaconda.sh && \
#     /bin/bash ~/anaconda.sh -b -p $CONDA_DIR && \
#     rm ~/anaconda.sh
# ENV PATH=$CONDA_DIR/bin:$PATH
# # 初始化 conda
# RUN conda init bash

# # 更新 conda 和所有包
# RUN conda update -n base -c defaults conda && \
#     conda update --all && conda clean -afy

ENV FASTAHACK $SOFT/fastahack/fastahack
# RUN wget https://anaconda.org/bioconda/fastahack/1.0.0/download/linux-64/fastahack-1.0.0-h4ac6f70_1.tar.bz2 && \
#     tar -jxvf fastahack-1.0.0-h4ac6f70_1.tar.bz2 && \
#     rm fastahack-1.0.0-h4ac6f70_1.tar.bz2 && \
#     mkdir -p $SOFT/fastahack && mv bin/fastahack $SOFT/fastahack/fastahack && \
#     rm -rf bin/ && rm -rf info/

RUN wget https://github.com/ekg/fastahack/archive/refs/tags/v1.0.0.tar.gz -O fastahack-1.0.0.tar.gz && \
tar -xzvf fastahack-1.0.0.tar.gz && \
rm fastahack-1.0.0.tar.gz && \
cd fastahack-1.0.0 && \
make && mkdir -p $SOFT/fastahack && mv fastahack $FASTAHACK && \
cd ../ && rm -rf fastahack-1.0.0/

ENV MS_V2 $SOFT/msisensor/db/microsatellites_selected.list
ADD microsatellites_selected.list $SOFT/msisensor/db/

ENV MSIFINDER $SOFT/MSIFinder-master/MSIdetect_v3.py
ADD *.py $SOFT/MSIFinder-master/
ADD requirements.txt $SOFT/MSIFinder-master/
RUN pip install -r $SOFT/MSIFinder-master/requirements.txt -i https://mirrors.aliyun.com/pypi/simple

RUN chmod +x $MSIFINDER

# 设置工作目录
WORKDIR /data

