import os
import csv
import argparse
from collections import defaultdict

def process_graph(input_file, output_dir, cri_threshold, target_city):
    # 確保輸出目錄存在
    if not os.path.exists(output_dir):
        os.makedirs(output_dir)

    # 打開原始 CSV 文件
    with open(input_file, mode='r', encoding='utf-8') as file:
        reader = csv.DictReader(file)

        graph_data = {}
        filtered_transactions = []

        for row in reader:
            buyer = row['buyer_masterid']
            seller = row['bidder_masterid']
            transaction_amount = row['tender_finalpriceUsd']
            cri_value = float(row['cri'])

            # 動態過濾條件：根據輸入參數判斷城市
            buyer_city = row['buyer_city']
            tender_year = int(float(row['tender_year']))
            tender_supplytype = row['tender_supplytype']
            tender_lotscount = row['tender_lotscount']

            # 如果 buyer_city 與 target_city 不匹配，跳過該條數據
            if buyer_city != target_city:
                continue
            if not (2015 <= tender_year <= 2019):
                continue
            if tender_supplytype not in ['WORKS', 'SERVICES']:
                continue
            if not tender_lotscount or int(tender_lotscount) == 0:
                continue
            if not buyer or not seller or not transaction_amount:
                continue

            transaction_count = int(tender_lotscount)
            filtered_transactions.append((buyer, seller, transaction_amount, cri_value, transaction_count))

        for index, (buyer, seller, transaction_amount, cri_value, transaction_count) in enumerate(filtered_transactions):
            if buyer not in graph_data:
                graph_data[buyer] = {'nodes': set(), 'edges': [], 'cri': cri_value}
            else:
                graph_data[buyer]['cri'] = max(graph_data[buyer]['cri'], cri_value)

            graph_data[buyer]['nodes'].add(buyer)
            graph_data[buyer]['nodes'].add(seller)

            if transaction_count == 1:
                label = 0
            elif 2 <= transaction_count <= 5:
                label = 1
            else:
                label = 2

            graph_data[buyer]['edges'].append((buyer, seller, label))

        # 生成圖文件
        graph_filename = f"{output_dir}/{dataset_prefix}_graph.txt" if cri_threshold == 0.5 else f"{output_dir}/{dataset_prefix}_cri{int(cri_threshold * 10)}_graph.txt"
        with open(graph_filename, mode='w') as file:
            for graph_id, (buyer, data) in enumerate(graph_data.items()):
                file.write(f"t # {graph_id}\n")
                node_id_map = {}
                for node_id, node in enumerate(data['nodes']):
                    node_id_map[node] = node_id
                    file.write(f"v {node_id} {node_id}\n")
                for (source, target, label) in data['edges']:
                    file.write(f"e {node_id_map[source]} {node_id_map[target]} {label}\n")

        # 生成標籤文件
        label_filename = f"{output_dir}/{dataset_prefix}_label.txt" if cri_threshold == 0.5 else f"{output_dir}/{dataset_prefix}_cri{int(cri_threshold * 10)}_label.txt"
        with open(label_filename, mode='w') as file:
            for buyer, data in graph_data.items():
                label = 1 if data['cri'] > cri_threshold else 0
                file.write(f"{label}\n")

        print(f"Graph saved to: {graph_filename}")
        print(f"Labels saved to: {label_filename}")

if __name__ == "__main__":
    parser = argparse.ArgumentParser(description="Process graph from CSV")
    parser.add_argument("-d", "--dataset", required=True, help="Dataset prefix, e.g., 'GE' or 'UK'")
    parser.add_argument("--cri", type=float, default=0.5, help="CRI threshold value, default is 0.5")
    args = parser.parse_args()

    dataset_prefix = args.dataset
    cri_threshold = args.cri

    # 動態設置 target_city
    target_city_mapping = {
        "GE": "Berlin",
        "UK": "London"
    }
    target_city = target_city_mapping.get(dataset_prefix, "Unknown")
    if target_city == "Unknown":
        print(f"Error: Unsupported dataset prefix '{dataset_prefix}'.")
        exit(1)

    input_file = f"{dataset_prefix}_DIB_2023.csv"
    print(input_file)

    # 根據 cri 決定資料夾名稱
    if cri_threshold == 0.5:
        output_dir = f"Pang/data/{dataset_prefix}"
    else:
        output_dir = f"Pang/data/{dataset_prefix}_cri{int(cri_threshold * 10)}"

    process_graph(input_file, output_dir, cri_threshold, target_city)
