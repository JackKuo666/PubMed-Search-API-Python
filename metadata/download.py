import json
import os
import ftplib
import hashlib
import logging
from time import sleep

# 配置参数
ftp_url = "download.nmdc.cn"
ftp_dir = "/pubmed"
local_dir = "./pubmed_data"
file_patterns = (".xml.gz", ".md5")  # 需要下载的文件类型
max_retries = 3  # 单个文件下载失败重试次数
retry_delay = 10  # 重试间隔(秒)

# 初始化日志
logging.basicConfig(
    filename='pubmed_download.log',
    level=logging.INFO,
    format='%(asctime)s - %(levelname)s - %(message)s',
    encoding='utf-8'
)

def verify_md5(file_path, expected_md5):
    """验证文件的MD5校验码"""
    hash_md5 = hashlib.md5()
    with open(file_path, "rb") as f:
        for chunk in iter(lambda: f.read(4096), b""):
            hash_md5.update(chunk)
    return hash_md5.hexdigest() == expected_md5

def download_file(ftp, filename, local_path):
    """带重试机制的文件下载"""
    for attempt in range(max_retries):
        try:
            with open(local_path, 'wb') as f:
                ftp.retrbinary(f"RETR {filename}", f.write)
            return True
        except Exception as e:
            logging.error(f"Attempt {attempt+1} failed for {filename}: {str(e)}")
            if attempt < max_retries - 1:
                sleep(retry_delay)
    os.remove(local_path)
    return False

def main():
    # 创建本地目录
    os.makedirs(local_dir, exist_ok=True)
    error_files = []
    try:
        # 连接FTP
        ftp = ftplib.FTP(ftp_url, timeout=30)
        ftp.login()  # 匿名登录
        ftp.cwd(ftp_dir)
        
        # 获取文件列表
        files = [f for f in ftp.nlst() if f.endswith(file_patterns) and f.startswith('pubmed25')]
        logging.info(f"Found {len(files)} files to process")

        # 下载主循环
        for idx, filename in enumerate(files, 1):
            try:
                local_path = os.path.join(local_dir, filename)
                
                # 跳过已存在的完整文件
                if os.path.exists(local_path):
                    logging.info(f"Skipping existing file: {filename}")
                    continue

                # 下载文件
                logging.info(f"Downloading ({idx}/{len(files)}) {filename}")
                success = download_file(ftp, filename, local_path)
                
                # MD5校验
                if success and filename.endswith(".md5"):
                    continue  # 不需要校验MD5文件本身
                else:
                    logging.info(f"Downloading MD5 for {filename}")
                    md5_file = filename + ".md5"
                    md5_local_path = os.path.join(local_dir, md5_file)
                    success2 = download_file(ftp, md5_file, md5_local_path)
                    
                if success and success2 and filename.endswith(".xml.gz"):
                    md5_file = filename + ".md5"
                    if md5_file in files:
                        md5_path = os.path.join(local_dir, md5_file)
                        with open(md5_path) as f:
                            expected_md5 = f.read().split()[1].replace('\n', '')
                        if not verify_md5(local_path, expected_md5):
                            logging.error(f"MD5 mismatch for {filename}")
                            os.remove(local_path)  # 删除损坏文件
            except Exception as e:
                logging.error(f"{filename}:Fatal error: {str(e)}")
                error_files.append(filename)

    except Exception as e:
        logging.error(f"Fatal error: {str(e)}")
    finally:
        if 'ftp' in locals():
            ftp.quit()
        with open('./error_file.json', 'w') as f:
            json.dump(error_files, f)

if __name__ == "__main__":
    main()
