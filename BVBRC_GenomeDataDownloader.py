# %%
# 導入會用到的模組
import csv, os, re, subprocess

# %%
# 定義一個變數，儲存 CSV 檔案的路徑
print("Please enter csv_file ABSOLUTE path: ")
csv_file = input()
print(csv_file, type(csv_file))

# 使用 csv 模組，讀取 CSV 檔案的第一行，並將其分割為四個欄位
with open(csv_file) as f:
  rows = csv.reader(f, delimiter=",")
  diseases = next(rows) # 取得第一行，即病症的名稱
  print("Groups ready to download:", diseases)

# 定義一個變數，儲存要下載檔案的 FTP 網址
ftp_url = "ftp://ftp.bvbrc.org/genomes/"

# 定義一個列表，儲存要下載的檔案格式
formats = [".PATRIC.gff", ".fna"]

# 顯示當前工作目錄
print("Current working directory:", os.getcwd())

# %%
# 對於每一種病症，建立一個資料夾，並進入
for i, disease in enumerate(diseases):

    os.mkdir(disease)
    os.chdir(disease)

    ## (1) 紀錄 disease 的 genome IDs 到所屬的資料夾
    with open(csv_file) as f, open(disease + "-genome-ids.txt", "w+") as ids:
        rows = csv.reader(f, delimiter=",")

        next(rows) # 跳過第一行，即 group names

        for row in rows:
            if row[i].strip(): # 如果 genome ID 不是空的
                ids.write(row[i].strip() + "\n") # 寫入 genome IDs

        # 移動 filestream 的指標到檔案開頭，重新讀檔
        ids.seek(0)
        rstrip_ids = ids.readlines()
        rstrip_ids[-1] = rstrip_ids[-1].rstrip() # 移除最後一行的換行符號

        ids.seek(0)
        ids.truncate()
        ids.writelines(rstrip_ids)
    ################################################################################################################
        

    ## (2) 對於指定的 format (gff 和 fna) 下載對應格式的檔案，並放進對應的資料夾
    for format in formats:
        
        # 建立存放對應格式檔案的資料夾，並進入
        os.mkdir(format)
        os.chdir(format)

        # 對於所有 genome ID，下載對應格式的檔案
        with open("../" + disease + "-genome-ids.txt") as ids:
            for id in ids:
                
                result = subprocess.run(["wget", ftp_url + id.strip() + "/" + id.strip() + format])

                if result.returncode != 0:
                    with open("error-ids.txt", "a") as error_ids:
                        error_ids.write(id.strip() + "\n")

                # os.system("wget " + ftp_url + id.strip() + "/" + id.strip() + format)

        # 離開 format 的資料夾，返回上一層
        os.chdir("..")
    ###################################################################################################################
        

    """
    若不需要這個功能，可以將 (3)、(4) 這兩段程式碼註解掉。
    
    以下的程式碼，使的最後產生的 GFF 檔案會包含 .fna 的內容；可以輸入 Roary 進行 pan-genome 分析。

    1) 會將 .fna 每條 contigs 的 FASTA 標題開頭加上 "accn|" 。
    (如：>accn|SPUA01000001   SPUA01000001.1   [Klebsiella pneumoniae strain NICU_1_P7 | 573.31029])
    
    2) 將所有 .gff 檔案從 "./.PATRIC.gff" 資料夾移動到上一層目錄，也就是當前的 disease 資料夾。
    
    3) 然後再將處理後的 .fna 內容附加到 .gff 文件最後。

    4) 最後會刪除 .fna 和 .PATRIC.gff 資料夾。
    
    這樣做是為了讓 .fna 和 .gff 檔案的內容能夠對應上。(.gff 的每個 feature id 都有一樣的前綴)
    
    """


    ## (3) 確認並取得成功下載 .gff 和 .fna 檔案的 genome IDs
    os.chdir("./.PATRIC.gff")

    # 取得在當前目錄下執行 BASH 指令 "ls" 的 standard output，並將其分割成一個 list
    files = subprocess.run(["ls"], capture_output=True, text=True)
    files_list = files.stdout.strip().split("\n")

    # 定義一個 regex，用來篩選特定格式的檔名
    regex_gff = re.compile(r"^\d+\.\d+\.PATRIC\.gff$")

    # 執行篩選
    gffs = []
    for file in files_list:
        if regex_gff.match(file):
            gffs.append(file)

    # 提取出檔名中的 genome ID
    gff_ids = []
    for name in gffs:
        ext_index = name.rfind(".PATRIC.gff")
        gff_id = name[:ext_index]
        gff_ids.append(gff_id)

    # 進到 ../.fna 資料夾，執行類似的操作
    os.chdir("../.fna")

    files = subprocess.run(["ls"], capture_output=True, text=True)
    files_list = files.stdout.strip().split("\n")

    regex_fna = re.compile(r"^\d+\.\d+\.fna$")

    fnas = []
    for file in files_list:
        if regex_fna.match(file):
            fnas.append(file)

    fna_ids = []
    for name in fnas:
        ext_index = name.rfind(".fna")
        fna_id = name[:ext_index]
        fna_ids.append(fna_id)

    # 返回上一層目錄: disease 資料夾
    os.chdir("..")

    # 比對 .gff 和 .fna 檔案的 genome IDs，找出兩者共有的 genome IDs；
    # 這些 IDs 將會用來一一執行下方合併檔案的操作
    ids_ready_to_merge = list(set(gff_ids) & set(fna_ids))
    ############################################################################################################################


    # (4) 合併 .fna 和 .gff 檔案
    for id in ids_ready_to_merge:

        # (4-1)  .fna 檔案中的每條 contigs，其 FASTA 標題將會從 ">" 被取代成 ">accn|"，並去除空格後的內容
        # eg. ">SPUA01000001   SPUA01000001.1   [Klebsiella pneumoniae strain NICU_1_P7 | 573.31029]" => ">accn|SPUA01000001"
        with open("./.fna/" + id.strip() + ".fna", "r+") as fna:

            # 將 fna 讀取為一個 list，每行為一個元素
            lines = fna.readlines()

            # 定義一個 generator，用來找出 FASTA 標題的 index
            def get_fasta_heading_index(lines):
                for i, line in enumerate(lines):
                    if line.startswith(">"):
                        yield i

            # 對於 generator 吐出的每一個 index，去除行末換行符號、只保留空白前部分、將 ">" 取代為 ">accn|"，最後再加上換行符號；
            # 這些更動暫時保存在 lines 這個 list 中
            for i in get_fasta_heading_index(lines):
                lines[i] = lines[i].rstrip().split(" ", 1)[0].replace(">", ">accn|") + "\n"
            
            # filestream 的指標移到檔案開頭，清空檔案內容；再將 lines 寫入檔案
            fna.seek(0)
            fna.truncate()
            fna.writelines(lines)
        # ----------------------------------------------------------------------------------------
        
        
        # (4-2) 找到對應的 .gff 檔案，append 處理後的 .fna 內容
        with open("./.PATRIC.gff/" + id.strip() + ".PATRIC.gff", "a") as gff:

            # 新增一行 "##FASTA" 標籤
            gff.write("##FASTA\n")

            # 將 .fna 的內容附加到 .gff 文件最後
            with open("./.fna/" + id.strip() + ".fna") as fna:
                gff.write(fna.read())
        
        # 移動 .gff 檔案到當前的 disease 資料夾
        subprocess.run(["mv", "./.PATRIC.gff/" + id.strip() + ".PATRIC.gff", "."])
        # ----------------------------------------------------------------------------------------

    # # 刪除 .fna 資料夾
    # os.system("rm -r ./.fna")
    
    # # 刪除 .PATRIC.gff 資料夾
    # os.system("rm -r ./.PATRIC.gff")                     
            
    ########################################################################################################################

    # 離開目前 disease 的資料夾，返回上一層，準備處理下一個 disease
    os.chdir("..")


