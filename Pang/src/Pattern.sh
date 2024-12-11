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
outputCGSPANSupport="../data/${dataset}/${dataset}_CGSupport.txt"
outputTKG="../data/${dataset}/${dataset}_TKG.txt"

# 使用 /mnt/c/... 路徑格式，此處要自己調整java.exe的路徑，可以把生成的txt檔案貼到_patter.txt底下進行操作
"/mnt/c/Users/pinyen/anaconda3/envs/test_env/Library/bin/java.exe" -Xmx8192m -jar spmf.jar run GSPAN "$input" "$outputCGSPAN" 0.9 8 true false true 
"/mnt/c/Users/pinyen/anaconda3/envs/test_env/Library/bin/java.exe" -Xmx8192m -jar spmf.jar run CGSPANSupport "$input" "$outputCGSPANSupport" 0.9 8 true false true
"/mnt/c/Users/pinyen/anaconda3/envs/test_env/Library/bin/java.exe" -Xmx8192m -jar spmf.jar run TKG "$input" "$outputTKG" 4 8 true false true
