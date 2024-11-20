#!/bin/bash

# 默認值
dataset=""

# 解析參數
while getopts "d:" opt; do
  case $opt in
    d) dataset="$OPTARG" ;;  # 將 -d 的值存入變量
    *) echo "Usage: $0 -d <dataset>"; exit 1 ;;  # 顯示用法並退出
  esac
done

# 檢查是否提供了 -d 參數
if [ -z "$dataset" ]; then
  echo "Error: -d parameter is required."
  exit 1
fi

# 設置輸入和輸出路徑
input="../data/${dataset}/${dataset}_graph.txt"
outputGSPAN="../data/${dataset}/${dataset}_pattern.txt"
outputCGSPAN="../data/${dataset}/${dataset}_CG.txt"

# 使用 /mnt/c/... 路徑格式
"/mnt/c/Users/Lockdream/miniconda3/envs/myenv3715/Library/bin/java.exe" -Xmx8192m -jar spmf.jar run GSPAN "$input" "$outputGSPAN" 0.9 8 true false true 
# "/mnt/c/Users/Lockdream/miniconda3/envs/myenv3715/Library/bin/java.exe" -Xmx8192m -jar spmf.jar run CGSPAN "$input" "$outputCGSPAN" 0.9 8 true false true