FROM registry.cn-shanghai.aliyuncs.com/kube-shipper/bash:5.1.4
RUN apk update && apk add --no-cache curl
RUN curl https://dl.min.io/client/mc/release/linux-amd64/mc -o /usr/bin/mc  &&  chmod +x /usr/bin/mc