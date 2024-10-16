# 使用带有Java环境的基础镜像
FROM openjdk:8-jdk
WORKDIR /data1
ENV SOFT /data1/workdir/software
ENV SAMTOOLS_1_8=$SOFT/samtools-1.8/samtools
ENV PICARD $SOFT/picard_v2.20.0/picard.jar
ENV GENCORE=$SOFT/gencore/gencore  

RUN apt-get update && apt-get install -y build-essential gcc g++ gdb clang make autoconf libtool libbz2-dev liblzma-dev libncurses5-dev libssl-dev zlib1g-dev  \
    && apt-get clean \
    && rm -rf /var/lib/apt/lists/*

RUN wget https://github.com/samtools/samtools/releases/download/1.8/samtools-1.8.tar.bz2 && \
    tar -xjvf samtools-1.8.tar.bz2 && \
    rm samtools-1.8.tar.bz2 && \
    cd samtools-1.8 && \
    ./configure --prefix=$SOFT/samtools-1.8 && \
    make && make install && \
    cd ../ && rm -rf samtools-1.8/

RUN wget http://sourceforge.net/projects/bio-bwa/files/bwa-0.7.12.tar.bz2 && \
    tar -jxvf bwa-0.7.12.tar.bz2 && \
    rm bwa-0.7.12.tar.bz2 && \
    cd bwa-0.7.12 && \
    make && \
    mv bwa /usr/local/bin/bwa && \
    cd ../ && rm -rf bwa-0.7.12/

RUN mkdir -p $SOFT/gencore && wget http://opengene.org/gencore/gencore -O $GENCORE && chmod a+x $GENCORE
  
RUN mkdir -p $SOFT/picard_v2.20.0 && wget https://github.com/broadinstitute/picard/releases/download/2.20.0/picard.jar -O $PICARD  &&  chmod +x $PICARD 
##&& java -jar $PICARD -h  

ENV SAMTOOLS_1_9=$SOFT/samtools-1.9/samtools
RUN wget https://github.com/samtools/samtools/archive/refs/tags/0.1.19.tar.gz -O samtools-0.1.19.tar.gz && \
    tar -xzvf samtools-0.1.19.tar.gz && \
    rm samtools-0.1.19.tar.gz && \
    cd samtools-0.1.19 && \
    make && mkdir -p $SOFT/samtools-1.9 && mv samtools $SAMTOOLS_1_9 && \
    cd ../ && rm -rf samtools-0.1.19/  

ENV PATH=$SOFT/samtools-1.8/bin:$PATH  

ADD aln_ctDNA.sh aln_post_ctDNA.sh /usr/local/bin/

# 赋予执行权限
RUN chmod +x /usr/local/bin/aln_ctDNA.sh && chmod +x /usr/local/bin/aln_post_ctDNA.sh

# 设置工作目录
WORKDIR /data

