# bio-info

### 输入文件

存放生信分析流程输入文件：[https://gitee.com/lingbohome/bio-data.git](https://gitee.com/lingbohome/bio-data.git)

### 容器内 Debug
```bash
docker run -it --rm -v /root/lingbo/bio/poc/data:/data  --name bio-demo --entrypoint bash registry.cn-shanghai.aliyuncs.com/bio-cloud/bioinfo-aln:v0.2
```

### 生成python依赖
```bash 
pipreqs .\06_MSI\ --force
```
