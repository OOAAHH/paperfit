#!/bin/bash

#判断是否有输入GSE ID号
if [ -z "$1" ]
then
        echo "Please provide a GSE number."
        exit 1
fi

#在当前目录下创建一个文件夹，以GSE ID号命名
mkdir -p $1
prefix=$(echo $1 | cut -c 1-$((${#1}-3)))

#递归下载该GSE下的所有元数据
proxychains4 wget -r -c -nH --cut-dirs=3 ftp://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}nnn/$1/ -P /Users/sunhao/DRaw/
#wget -r -c -nH --cut-dirs=3 ftp://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}nnn/$1/ -P /Users/sunhao/DRaw/

#输出下载成功
echo "Successed downloading $1"
