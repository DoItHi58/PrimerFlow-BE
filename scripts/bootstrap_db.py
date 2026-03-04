import os
import sys
import urllib.parse
import urllib.request

CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)

DEFAULT_DB_PATH = os.path.join(BASE_DIR, "database", "annotations.db")
DEFAULT_DB_DIR = os.path.dirname(DEFAULT_DB_PATH)
DEFAULT_GENOME_PATH = os.path.join(BASE_DIR, "database", "raw_data", "GRCh38.primary_assembly.genome.fa.gz")
DOWNLOAD_CHUNK_SIZE = 1024 * 1024

def _env(name, default=None):
    value = os.getenv(name)
    if value:
        return value
    return default

def _ensure_dir(path):
    if path:
        os.makedirs(path, exist_ok=True)

def _download_stream(resp, dest_path):
    with open(dest_path, "wb") as out_file:
        while True:
            chunk = resp.read(DOWNLOAD_CHUNK_SIZE)
            if not chunk:
                break
            out_file.write(chunk)

def _download_file(url, dest_path):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp:
        _download_stream(resp, dest_path)

def _download_google_drive(url, dest_path):
    req = urllib.request.Request(url, headers={"User-Agent": "Mozilla/5.0"})
    with urllib.request.urlopen(req) as resp:
        cookies = resp.headers.get_all("Set-Cookie") or []
        confirm_token = None
        for cookie in cookies:
            if "download_warning" in cookie:
                parts = cookie.split(";")[0].split("=")
                if len(parts) == 2:
                    confirm_token = parts[1]
                    break

        if not confirm_token:
            _download_stream(resp, dest_path)
            return

    parsed = urllib.parse.urlparse(url)
    query = urllib.parse.parse_qs(parsed.query)
    query["confirm"] = [confirm_token]
    new_query = urllib.parse.urlencode(query, doseq=True)
    url = urllib.parse.urlunparse(parsed._replace(query=new_query))
    _download_file(url, dest_path)


def _download_any(url, dest_path):
    if "drive.google.com" in url:
        _download_google_drive(url, dest_path)
    else:
        _download_file(url, dest_path)


def _download_if_missing(url, dest_path, label):
    if os.path.exists(dest_path):
        print(f"{label} 파일이 이미 존재합니다. 스킵: {dest_path}")
        return
    if not url:
        return
    _ensure_dir(os.path.dirname(dest_path))
    print(f"{label} 다운로드 시작: {url}")
    _download_any(url, dest_path)
    if not os.path.exists(dest_path):
        print(f"{label} 다운로드 실패: 파일이 생성되지 않았습니다.")
        sys.exit(1)
    print(f"{label} 다운로드 완료: {dest_path}")


def _validate_bgzip_sidecars(genome_path):
    if not genome_path.endswith(".gz"):
        return
    fai_path = f"{genome_path}.fai"
    gzi_path = f"{genome_path}.gzi"
    missing = [path for path in (fai_path, gzi_path) if not os.path.exists(path)]
    if missing:
        print("FASTA 인덱스 파일이 누락되었습니다.")
        print("- GENOME_PATH가 .gz일 때는 bgzip FASTA + .fai + .gzi가 모두 필요합니다.")
        print(f"- 누락 파일: {', '.join(missing)}")
        print("- 해결: GENOME_FAI_URL / GENOME_GZI_URL 환경변수를 설정하거나 파일을 사전에 배치하세요.")
        sys.exit(1)

def main():
    db_url = _env("DB_URL")
    db_path = _env("DB_PATH", DEFAULT_DB_PATH)
    db_dir = os.path.dirname(db_path)
    genome_url = _env("GENOME_URL")
    genome_path = _env("GENOME_PATH", DEFAULT_GENOME_PATH)
    genome_fai_url = _env("GENOME_FAI_URL")
    genome_gzi_url = _env("GENOME_GZI_URL")
    genome_dir = os.path.dirname(genome_path)

    if not db_url:
        print("DB_URL 환경변수가 없습니다.")
        sys.exit(1)

    if os.path.exists(db_path):
        print(f"DB 파일이 이미 존재합니다. 스킵: {db_path}")
    else:
        _ensure_dir(db_dir)
        print(f"DB 다운로드 시작: {db_url}")
        _download_any(db_url, db_path)
        if not os.path.exists(db_path):
            print("DB 다운로드 실패: 파일이 생성되지 않았습니다.")
            sys.exit(1)
        print(f"DB 다운로드 완료: {db_path}")

    if genome_url:
        if os.path.exists(genome_path):
            print(f"FASTA 파일이 이미 존재합니다. 스킵: {genome_path}")
        else:
            _ensure_dir(genome_dir)
            print(f"FASTA 다운로드 시작: {genome_url}")
            _download_any(genome_url, genome_path)
            if not os.path.exists(genome_path):
                print("FASTA 다운로드 실패: 파일이 생성되지 않았습니다.")
                sys.exit(1)
            print(f"FASTA 다운로드 완료: {genome_path}")
    else:
        print("GENOME_URL 환경변수가 없습니다. FASTA 다운로드를 건너뜁니다.")

    if os.path.exists(genome_path) and genome_path.endswith(".gz"):
        _download_if_missing(genome_fai_url, f"{genome_path}.fai", "FAI")
        _download_if_missing(genome_gzi_url, f"{genome_path}.gzi", "GZI")
        _validate_bgzip_sidecars(genome_path)

if __name__ == "__main__":
    main()
