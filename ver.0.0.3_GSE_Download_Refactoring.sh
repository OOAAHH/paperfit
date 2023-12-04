#!/bin/bash
# 需要提前部署 proxychains, parallel 工具, 并对proxychains.conf进行设置

# 关于for循环中的使用变量的写法是否有花括号
# 由于变量后面是一个斜杠 /，这个斜杠已经清楚地标示了变量名的结束，
# 所以花括号是可选的。因此，"$CHECK_DIR"/* 和 "${CHECK_DIR}"/* 在功能上是相同的，都是为了获取CHECK_DIR目录下的所有文件。

# 全局变量 请勿更改
DOWNLOAD_URL="ftp://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}nnn/$1/suppl/"
CHECK_URL="https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}nnn/$1/suppl/"
DOWNLOAD_DIR="$2"
CHECK_DIR="$2/$1/suppl"
LOG_DIR="$2/proxychains_download_log"
# 全局变量 用户自定义log文件
INCONSISTENT_FILES="$2/inconsistent_files.txt"
SIZE_DIFFERENCE_FILE="$2/diff.txt"
LOG_FILE="$2/log.txt"
ERROR_LOG="$2/download_errors.log"
# 全局变量 用户自定义并行下载数
MAX_JOBS=5

# 创建目录 可有可无, 按照此处使用的wget下载命令、以及NCBI ftp服务器的命名规则, 
# 以GSE编号为文件名的文件夹会自动被创建
create_directory() {
    mkdir -p "$1"
}


download_file() {
    local gse_id=$1
    local prefix=$(echo $gse_id | cut -c 1-$((${#gse_id}-3)))
    local download_url="ftp://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}nnn/${gse_id}"
    proxychains4 wget -r -c -nH --cut-dirs=3 "$download_url" -P "$DOWNLOAD_DIR" -o "$LOG_FILE" 2>&1 &
    # show_spinner $!
    if [ $? -ne 0 ]; then
        echo "wget下载'$gse_id'失败，重试次数超标。" >> "$ERROR_LOG"
    else
        echo "下载完成"
        for local_file in "$CHECK_DIR"/*; do
            echo "开始校验 $local_file "
            check_file_size "$1"
        done
    fi
}

parallel_download() {
    local gse_list_file=$1
    local download_path=$2
    local max_jobs=$MAX_JOBS  # 从函数参数获取并行作业的数量
    # 检查输入文件是否存在
    if [ ! -f "$gse_list_file" ]; then
        echo "GSE ID list file not found: $gse_list_file"
        exit 1
    fi
    # 使用GNU parallel进行并行下载
    # 这里我们直接将GSE ID作为参数传递给脚本 $0 是调用函数自身
    cat "$gse_list_file" | parallel -j "$max_jobs" bash "$0" {} "$download_path"
}

redownload_file() {
    local filename=$1
    # 删除原有文件并尝试重新下载
    rm "$CHECK_DIR/$filename"
    proxychains4 wget -c "$DOWNLOAD_URL$filename" -P "$CHECK_DIR"
    # 检查下载是否成功
    if [ $? -eq 0 ]; then
        local new_size=$(stat -c%s "$local_dir/$filename")
        if [ "$ftp_size" != "$new_size" ]; then
            echo "$filename" >> $INCONSISTENT_FILES
            local difference=$(echo "scale=2; 100 * ($ftp_size - $new_size) / $ftp_size" | bc)
            echo "文件'$filename'大小差异百分比: $difference%"
            echo "文件'$filename'大小差异百分比: $difference%" >> $SIZE_DIFFERENCE_FILE
            echo "文件'$filename'已重新下载"
        fi
    else
        echo "文件'$filename'重新下载失败" >> $INCONSISTENT_FILES
    fi
}

check_file_size() {
    local gse_id=$1
    local prefix=$(echo $gse_id | cut -c 1-$((${#gse_id}-3)))
    local check_url="https://ftp.ncbi.nlm.nih.gov/geo/series/${prefix}nnn/$1/suppl/"
    for LOCAL_FILE in "$CHECK_DIR"/*; do
        #提取文件名,这里基于本地文件获取列表, 是我们假定了目录传输是完整的,只是个别的文件下载的不完整. 也可以改为从filelist文件进行循环
        FILENAME=$(basename "$LOCAL_FILE")
        SERVER_RESPONSE=$(wget --spider --server-response "$check_url$FILENAME" 2>&1)
        echo "服务器响应: $SERVER_RESPONSE"
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
            #删除该文件重新下载 使用函数redownload_file()
            redownload_file "$FILENAME"
        else
            echo "文件下载成功"
            for LOCAL_FILE in "$CHECK_DIR"/*; do
                FILENAME=$(basename "$LOCAL_FILE")
                if [[ "$LOCAL_FILE" =~ \.gz$ ]]; then
                    if ! gunzip -t "$LOCAL_FILE"; then
                        echo "压缩文件 '$LOCAL_FILE' 损坏，重新下载。" >> "$ERROR_LOG"
                        #删除该文件重新下载 使用函数redownload_file()
                        redownload_file "$FILENAME"
                    fi
                elif [[ "$LOCAL_FILE" =~ \.tar$ ]]; then
                    if ! tar -tf "$LOCAL_FILE" > /dev/null; then
                        echo "压缩文件 '$LOCAL_FILE' 损坏，重新下载。" >> "$ERROR_LOG"
                        #删除该文件重新下载 使用函数redownload_file()
                        redownload_file "$FILENAME"
                    fi
                fi
            done
        fi
    done
}

# 动态显示函数
show_spinner() {
    local pid=$1
    local delay=0.1
    local spinstr='|/-\'
    while [ "$(ps a | awk '{print $1}' | grep $pid)" ]; do
        local temp=${spinstr#?}
        printf " [%c]  " "$spinstr"
        local spinstr=$temp${spinstr%"$temp"}
        sleep $delay
        printf "\b\b\b\b\b\b"
    done
    printf "    \b\b\b\b"
}


# 主函数
main() {
    # 判断键盘输入值 $1 是GSE编号或者GSE编号的文件 $2 是下载路径
    if [ -z "$1" ] || [ -z "$2" ]; then
        echo "Usage: $0 <GSE number> <path>"
        exit 1
    fi
    # 如果输入的是列表, 则使用并行下载 这里特指只包含GSE编号的文件, 这个文件最好以单列呈现, 因为实际上调用了cat命令
    if [ -f "$1" ]; then
        parallel_download "$1" "$2" "$MAX_JOBS"
        # 如果不是就启动单一下载
    else
        download_file "$1" "$2"
    fi
}

main "$@"
