# read_file.py
'''
使用conda 环境 

GEO数据库的一些格式, 暂时没有考虑loom文件
文件格式: csv, python
文件格式: gz, shell
文件格式: tar, shell
文件格式: txt, python
文件格式: h5, python
文件格式: h5ad, python
文件格式: xlsx, python
文件格式: bw, 数量: python
文件格式: rds/RDS, R
文件格式: h5Seurat, R
文件格式: zip, python
文件格式: tsv, python
文件格式: mtx, python
'''

import sys
import pandas as pd
import gzip
import h5py
import pyBigWig
import scanpy
import os
os.environ['R_HOME'] = 'you_conda_path/lib/R' # from conda env R423
import rpy2.robjects as ro
from rpy2.robjects import pandas2ri


# for RDS
def read_rds(file_path):
    pandas2ri.activate()
    readRDS = ro.r['readRDS']
    return pandas2ri.ri2py(readRDS(file_path))
# for h5Seurat
def read_h5seurat(file_path):
    pandas2ri.activate()
    seurat = importr('Seurat')
    seurat_disk = importr('SeuratDisk')
    return seurat_disk.LoadH5Seurat(file_path)
# for others
def read_file(file_path, log_file):
    file_format = file_path.split('.')[-1].lower()
    try:
        if file_format in ['csv', 'tsv']:
            return pd.read_csv(file_path)
        elif file_format == 'txt':
            with open(file_path, 'r') as file:
                return file.read()
        elif file_format == 'xlsx':
            return pd.read_excel(file_path)
        elif file_format == 'gz':
            with gzip.open(file_path, 'rb') as file:
                return file.read()
        elif file_format in ['h5', 'h5ad']:
            with h5py.File(file_path, 'r') as file:
                return file
        elif file_format == 'bw':
            with pyBigWig.open(file_path) as bw:
                return bw
        elif file_format == 'mtx':
            return scanpy.read_mtx(file_path)
        elif file_format in ['rds', 'rds']:
            return read_rds(file_path)
        elif file_format == 'h5seurat':
            return read_h5seurat(file_path)
        else:
            with open(log_file, 'a') as log:
                log.write(f"不支持的文件格式: {file_format}\n")
            return None
    except Exception as e:
        with open(log_file, 'a') as log:
            log.write(f"读取文件时出错: {file_path}, 错误: {e}\n")
        return None

if __name__ == "__main__":
    if len(sys.argv) < 3:
        print("Usage: python read_file.py <file_path> <log_file>")
    else:
        read_file(sys.argv[1], sys.argv[2])
