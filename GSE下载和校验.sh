#!/bin/bash
#忽略了一点, 加速节点本身也是有带宽限制的, 并行下载数目需要稍微留下冗余, 便于程序在retry时能够链接服务器
#判断是否有输入GSE ID号
if [ -z "$1" ]
then
    echo "Please provide a GSE number."
    exit 1
fi

#在当前目录下创建一个文件夹，以GSE ID号命名
mkdir -p $1
prefix=$(echo $1 | cut -c 1-$((${#1}-3)))
#FTP服务器的URL,进行FTP按照目录下载, 必须使用“ftp://”, 而这里进行查询的时候必须使用“https://”, 否则进行查询时获取的信息是不正确的, 这很奇怪, 或许和传输协议的具体原理有关.服务器方面为两种协议设置的不同的反应.
FTP_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}nnn/$1/suppl/"

SIZE_DIFFERENCE_FILE=/YOUPATH
#本地文件存储的目录
LOCAL_DIR="/YOUPATH/$1/suppl/"
#存储不一致文件列表的文件
INCONSISTENT_FILES="/YOUPPATH/inconsistent_files.txt"

#清空不一致文件列表
#> $INCONSISTENT_FILES

#递归下载该GSE下的所有元数据
#proxychains4 wget -r -c -nH --cut-dirs=3 ftp://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}nnn/$1/ -P /YOUPATH/
#wget -r -c -nH --cut-dirs=3 ftp://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}nnn/$1/ -P /YOUPATH/

proxychains4 wget -r -c -nH --cut-dirs=3 ftp://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}nnn/$1/ -P /YOUPATH/ -o /YOUPATH/log.txt
#根据wget返回值进行判断, 不过高设置重试次数的原因: 不希望出现意外的死循环
if [ $? -ne 0 ]
then
    echo "wget下载'$1'失败，重试次数超标。" >> /YOUPATH/download_errors.log
fi

#遍历本地目录中的文件
for LOCAL_FILE in "$LOCAL_DIR"/*; do
    #提取文件名,这里基于本地文件获取列表, 是我们假定了目录传输是完整的,只是个别的文件下载的不完整. 也可以改为从filelist文件进行循环
    FILENAME=$(basename "$LOCAL_FILE")
    SERVER_RESPONSE=$(wget --spider --server-response "$FTP_URL$FILENAME" 2>&1)
    #echo "服务器响应: $SERVER_RESPONSE"
    FTP_SIZE=$(echo "$SERVER_RESPONSE" | grep -i 'Content-Length:' | sed -e 's/.*Content-Length: \([0-9]*\).*/\1/')
    echo "服务器'$FILENAME'文件大小为:$FTP_SIZE"
    #获取本地文件的大小
    LOCAL_SIZE=$(stat -c%s "$LOCAL_FILE")
    echo "本地'$FILENAME'文件大小:$LOCAL_SIZE"
    #比较本地文件大小和FTP上的文件大小
    if [ "$FTP_SIZE" != "$LOCAL_SIZE" ]; then
        #如果大小不一致，将文件名写入列表,这里不直接启动下载,是希望避免出现意外的死循环
        echo "$FILENAME" >> $INCONSISTENT_FILES
        DIFFERENCE=$(echo "scale=2; 100 * ($FTP_SIZE - $LOCAL_SIZE) / $FTP_SIZE" | bc)
        echo "文件'$FILENAME'大小差异百分比: $DIFFERENCE%"
        echo "文件'$FILENAME'大小差异百分比: $DIFFERENCE%" >> $SIZE_DIFFERENCE_FILE
        #删除该文件重新下载
        echo "文件'$FILENAME'重新下载"
        rm "$LOCAL_FILE"
        DOWNLOAD_URL="$FTP_URL$FILENAME"
        proxychains4 wget -c "$DOWNLOAD_URL" -P "$LOCAL_DIR"

        #检查下载是否成功
        if [ $? -eq 0 ]; then
            LOCAL_SIZE=$(stat -c%s "$LOCAL_FILE")
            if [ "$FTP_SIZE" != "$LOCAL_SIZE" ]; then
                echo "$FILENAME" >> $INCONSISTENT_FILES
            DIFFERENCE=$(echo "scale=2; 100 * ($FTP_SIZE - $LOCAL_SIZE) / $FTP_SIZE" | bc)
            echo "文件'$FILENAME'大小差异百分比: $DIFFERENCE%"
            echo "文件'$FILENAME'大小差异百分比: $DIFFERENCE%" >> $SIZE_DIFFERENCE_FILE
            echo "文件'$FILENAME'已重新下载"
            fi
        else
            echo "文件'$FILENAME'重新下载失败" >> $INCONSISTENT_FILES
        fi
    fi
done

#检查不一致的文件列表是否为空
if [ -s $INCONSISTENT_FILES ]; then
    echo "以下文件大小不一致:"
    cat "$INCONSISTENT_FILES"
else
    echo "所有文件均已校验，无不一致。"
fi

echo "Successed downloading $1"
