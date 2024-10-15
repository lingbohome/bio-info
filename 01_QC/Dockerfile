# 使用带有Java环境的基础镜像
FROM openjdk:8-jdk
# 安装Trimmomatic
RUN apt-get update && apt-get install -y wget && \
    wget http://www.usadellab.org/cms/uploads/supplementary/Trimmomatic/Trimmomatic-0.36.zip && \
    unzip Trimmomatic-0.36.zip && \
    rm Trimmomatic-0.36.zip && \
    mv Trimmomatic-0.36 /usr/local/bin
RUN curl https://dl.min.io/client/mc/release/linux-amd64/mc -o /usr/bin/mc  &&  chmod +x /usr/bin/mc    
# 设置环境变量
ENV TRIMMOMATIC /usr/local/bin/Trimmomatic-0.36/trimmomatic-0.36.jar
ENV ADAPTER /path/to/adapter/sequences.fasta

# 将脚本复制到容器中
COPY qc_ctDNA.sh /usr/local/bin/qc_ctDNA.sh

# 赋予执行权限
RUN chmod +x /usr/local/bin/qc_ctDNA.sh

# 设置工作目录
WORKDIR /data

