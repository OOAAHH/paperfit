paperfit
这是用来批量寻找论文信息的代码，可用通过DOI号、GSE编号、PMID来获取文献基本信息


2023年11月4日

增加了GSE下载的shell脚本

# 在centOS上实现并行下载:
  1. 首先安装EPEL仓库

  `sudo yum install epel-release`

  2. 安装'proxychains'
     
  `sudo yum install proxychains-ng`

  3. 编辑配置文件'proxychains.conf'

  `sudo vi /etc/proxychains.conf`
  ###### 这里需要注意, 需要根据使用的是http或者socks5来具体的设置代理, 以socks5为例, 我的代理软件设置在7890端口进行监听. 那么本地代理的设置则需要你在`proxychains.conf`的末尾添加这样一行代码.

  `socks5  127.0.0.1 7890`

  4. 测试

  `proxychains4 wget ftp://ftp.example.com`
  如下图所示, `正在连接 ftp.ncbi...`这一行显示了, 本机和远方服务器224.0.0.1:22的连接被`proxychains`严格按照我们设置的socks5监听端口`127.0.0.1:7890`进行转发.
  ﻿![图片](https://github.com/OOAAHH/paperfit/assets/19518905/33582384-5cc0-45d1-8400-0398a65f34a9)
######   btw 通过不同的加速节点进行下载时, 速度有较大的差异. 目前尝试了多个供应商的节点, iggscholar的US节点有最快的下载速度: 14.5MB/s, 但magicloud拥有最多的每月不限速流量: 100TB.
  
  5. 安装并行工具
GNU Parallel 是一个 shell 工具，用于在一台机器上并行执行命令.

  `sudo yum install parallel`

  6. 测试并行下载

  这里我预备了一个txt文件:
|------------------|
| `./1.sh GSE194122` |
| `./1.sh GSE140203` |
| `./1.sh GSE164378` |
| `./1.sh GSE232523` |
| `./1.sh GSE182475` |

  每一行都是可以直接执行的shell命令, 为了实现并行下载, 我们可以使用`cat shell.txt | parallel -j 20`. 这里的`-j 20`即代表并行20个下载进程的意思, 一般来说parallel的并行的设置是一个核对应一个进程, 但在下载任务中, 还需要根据你的网速来具体的考虑并行进程的数量, 以实现下载成功率和下载速度的最优配置.








