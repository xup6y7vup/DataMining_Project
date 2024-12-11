
# 專案文件

## 專案結構

```plaintext
專案根目錄
├── Pang
│   ├── data       # 專案使用的資料集
│   │   ├── DD     
│   │   ├── FOPPA 
│   │   ├── GE     # 德國的公共標案購買
│   │   ├── MUTAG
│   │   ├── NCl1
│   │   ├── NCl109
│   │   ├── PRO    # 網路上已經標註好的蛋白質資料集
│   │   ├── PTC
│   │   └── UK     # 英國的公共標案購買
│   │ 
│   ├── result  # 結果產生的資料夾
│   │
│   └── src     # 主要程式的資料夾
│       ├── count_graph.py          # 計算有幾個圖     
│       ├── ECML.py                 # 使用論文預先設定好的資料集，用於重現論文
│       ├── PANG.py                 # 使用外部資料集時使用，用於實驗
│       ├── Pattern.sh              # 用於產生XXX_Pattern.txt
│       ├── ProcessingPatterns.py   # 用XXX_graph.txt 以及 XXX_Pattern.txt 產生 XXX_Mono.txt 以及 XXX_iso.txt
│       └── spmf.jar # 提供不同的演算法，可以設定演算法
│
├── environment.yml
├── gen_graph_lable.py    # 用於產生 XXX_graph.txt 以及 XXX_label.txt
└── README.md             # 專案說明文件
```

## 使用說明

### 1. 安裝相依套件

#### 安裝套件

1. 執行以下指令來安裝相依套件：
    ```bash
    conda env create -f environment.yml
    ```
    **注意**：請確定已經有安裝conda，可以在終端輸入 `conda --version` 確認是否安裝完成。 `environment.yml` 已經指定 Python 版本，但請於安裝完後確認版本是否正確，若顯示 Python 3.7.1 即可執行程式，環境名稱預設是 `name: test_env`，如果需要修改可以自行去 `environment.yml` 進行修改。

2. 執行以下指令來啟動環境：
    ```bash
    conda activate test_env
    ```

## 執行程式

### 1. 重現論文

請先 `cd` 到 `data` 資料夾中，將`NCI1`以及`NCI109`解壓縮。

請先 `cd` 到 `src` 資料夾中，然後於終端執行下列指令。如果沒有問題，會在 `result` 資料夾中看到 `table2`、`table4`、`table5` 和 `table6`，分別對應論文中的表格。

```bash
python ECML.py
```

**注意**：在 `ECML.py` 的 857~859 行有三種不同的 `results`，分別代表不同的格式。

```python
results = CORK  # 預設設定，建議使用此模式
results = Baseline
results = DGCNN
```

- **`results = CORK`** 是論文中的方法 (GenBin, GenOcc, IndBin, IndOcc, CloBin, CloOcc) 加上 `CORK` 的結果。
- **`results = Baseline`** 是論文中的方法加上 WL、WLOA、G2V 的結果。
- **`results = DGCNN`** 是論文中的方法加上 DGCNN 的結果。

**建議**：如果僅重現論文結果，請使用 `results = CORK`（已預設）。

### 2. 重現實驗

請先 `cd` 到 `src` 資料夾中，然後於終端輸入執行下列指令。

```bash
python PANG.py -d XXX -k XX
```

- **`-d`** 表示資料集名稱，例如 `UK`。
- **`-k`** 是論文中可調整的參數。

例如：
```bash
python PANG.py -d UK -k 10
```
執行後會在 `result` 資料夾中產生 `UK_10results.txt`，可以看到實驗的準確度。

## 自製資料集說明

此部分以 `-d UK -cri 0.7` 為範例說明。

### 1. 準備資料集

- 資料集下載鏈接：[GE_DIB_2023.csv 和 UK_DIB_2023.csv (Google Drive)](https://drive.google.com/drive/u/0/folders/1mQ_Imk1bTMJYn-1zlPcYPHVZdm6qYEaf)
- **注意**：解壓縮後請將資料與 `gen_graph_lable.py` 放在同目錄下。

### 2. 產生 `graph` 和 `label` 文件

執行下列指令以產生 `UK_cri7_graph.txt` 和 `UK_cri7_label.txt`：
```bash
python gen_graph_lable.py -d UK --cri 0.7
```

**注意**：
- `cri` 預設是 `0.5`，如果沒有輸入 `cri`，會產生 `UK` 資料夾，並生成 `UK_graph.txt` 和 `UK_label.txt`。

### 3. 產生 `pattern` 文件

`cd` 進 `Pang/src` 資料夾，執行下列指令以產生 `UK_cri7_pattern.txt`：
```bash
bash Pattern.sh -d UK_cri7
```

**注意**：
- 在 `Pattern.sh` 中，需確認 Java 的路徑，預設是我的安裝路徑，麻煩修改成自己java的安裝路徑，例如：
    ```bash
    "/mnt/c/Users/Lockdream/miniconda3/envs/myenv3715/Library/bin/java.exe" -> "your/path/jave.exe"
    ```

### 4. 產生 `Mono` 和 `Iso` 文件

執行下列指令以產生 `UK_cri7_Iso.txt` 和 `UK_cri7_Mono.txt`：
```bash
python ProcessingPatterns.py -d UK_cri7
```

### 5. 產生實驗結果

執行下列指令以產生結果。如果沒有問題，會在 `result` 資料夾中產生 `UK_cri7_10results.txt`：
```bash
python PANG.py -d UK_cri7 -k 10
```
