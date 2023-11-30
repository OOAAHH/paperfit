## paperfit
这是用来批量寻找论文信息的代码，可用通过DOI号、GSE编号、PMID来获取文献基本信息


### 2023年11月4日

增加了GSE下载的shell脚本,以及进行文件核对的ipynb文件

### 2023年11月8日

学到了wget的具有校验文件完整性的内建功能, 我为下载脚本增加了文件校验的步骤, 现在下载是否完成已经可以被标记了.

### 2023年11月10日

更改了下载文件的判断, 有些文件下载失败是对方服务器的问题, 莫得办法, 只能先记下了, 然后过一段时间再下载. 


### 2023年11月30日 重磅更新!!!

本次更新整理了原来的代码逻辑, 用函数的形式重构了脚本, 现在这个脚本原生支持并行下载! 原生支持下载文件后的校验! 很多时候GEO数据库的压缩文件在传输过程中会出错, 所以需要在下载完成后校验一下, 这在你需要下载很多数据的时候非常有用!

原本也不会更新的, 毕竟这个月太忙了, 我的leader在最初讨论的时候, 认为我们在生成环境部署了下载加速工具之后, 下载数据这个事情就不需要专门的人来做了, 这件事会变成“谁用谁下载”. 
但俗话说
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
“复杂性不会消失, 它只会转移”,
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~
在这件事上也是一样的. 好吧, 让我们有始有终.

学了面向对象编程的基本思想, 一直没地方用, 把这段时间的一些想法在这个下载脚本上实现了, 月初给自己画的饼🫓终于吃上了一部分. 🤔等我整明白GPT的api了就把🫓补全.

新的代码看起来真的`像诗一样优美`, 我竟然能搞出来这种东西! 此处着重推荐《代码整洁之道》这本书, 非常的amazing!!!
`[代码整洁之道.pdf](https://github.com/OOAAHH/paperfit/files/13513016/default.pdf)` 
##### 显然, 我这里分享的这个PDF文件仅供交流学习使用, 请在下载后24小时之内删除.

<img width="744" alt="截屏2023-11-30 21 16 10" src="https://github.com/OOAAHH/paperfit/assets/19518905/a3e45674-fad2-439b-887b-7b1d5a107e83">

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

“ ———Bjarne显然认为整洁的代码读起来令人愉悦. 读这种代码, 就像见到手工精美的音乐盒或者设计精良的汽车一般, 让你会心一笑. ”

                                                    ——Bjarne Stroustrup, C++ Programming Language 作者

~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~


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


# 进行文件核对
  1. 一般来说, 我们验证文件是否完整的思路是: 首先检查文件是否齐全, 再校验本地文件和服务器文件的md5值是否一致,
  而在从NCBI的GEO数据库的ftp服务器进行下载时, 很多情况下没有给出md5值,这时我的文件校验的思路是, 首先获取本地文件列表, 再通过python bs4爬虫从ftp服务器获取对应的文件结构和信息.
  进行文件比对的桥接信息是GSE编号, 本地的文件和ftp的文件均按照GSE编号进行存储. 我们可以获得如下图所示的文件信息列表, 此时我们就有了进行比对的基础, 使用代码或者wps/office均可, 我在这里展示的是wps的可视化结果. ps: 我小学的时候Excel的最大行数就是65536, 几十年了还是这样(bushi.
  <img width="1002" alt="截屏2023-11-04 17 14 11" src="https://github.com/OOAAHH/paperfit/assets/19518905/4c8463d6-63cf-48cc-9894-f2194b231c41">


  2. 另一个方法
  我们使用wget的`-c`参数来实现文件的校验, 它的功能是进行断点续传. 只要我们“重新下载”一次即可再一定程度上实现文件校验的目的. 它会遍历每个文件, 并不算是“聪明”的做法, 但是有效.
  <img width="311" alt="截屏2023-11-04 17 22 41" src="https://github.com/OOAAHH/paperfit/assets/19518905/ef373ba5-4675-486e-bae0-2ef6e3f273e3">


  ###### 2023年11月8日
  
  我本意是这个东西想想就得了, 所以简单的做了点东西, 想来我们的新下载节点不可能出这种问题的.
  BUT, 即使retry设置了20次, 还是有很多文件下载不完全, 离谱. 看日志发现可能是并行太多, 有些retry根本没有足够的带宽再连接服务器.


