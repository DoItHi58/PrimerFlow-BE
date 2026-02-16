import gzip
import os
import sqlite3

# ---------------------------------------------------------
# 1. 경로 및 설정
# ---------------------------------------------------------
CURRENT_DIR = os.path.dirname(os.path.abspath(__file__))
BASE_DIR = os.path.dirname(CURRENT_DIR)

DB_PATH = os.path.join(BASE_DIR, "database", "annotations.db")
RAW_DATA_DIR = os.path.join(BASE_DIR, "database", "raw_data")

# 제한 효소 목록 (필요시 추가)
ENZYMES = {
    'EcoRI': 'GAATTC',
    'BamHI': 'GGATCC',
    'HindIII': 'AAGCTT',
    'NotI': 'GCGGCCGC'
}

def get_db_connection():
    return sqlite3.connect(DB_PATH)

def init_schema(conn):
    cursor = conn.cursor()
    print("⚙️  스키마 재설정 중 (기존 데이터 삭제)...")
    
    # 기존 테이블 삭제 후 재생성 (Clean Build)
    cursor.executescript("""
        DROP TABLE IF EXISTS snp;
        DROP TABLE IF EXISTS restriction_site;
        DROP TABLE IF EXISTS exon;
        DROP TABLE IF EXISTS repeats;

        CREATE TABLE snp (id INTEGER PRIMARY KEY, chrom TEXT, pos INTEGER);
        CREATE TABLE restriction_site (id INTEGER PRIMARY KEY, chrom TEXT, name TEXT, start INTEGER, end INTEGER);
        CREATE TABLE exon (id INTEGER PRIMARY KEY, chrom TEXT, start INTEGER, end INTEGER, transcript_id TEXT);
        CREATE TABLE repeats (id INTEGER PRIMARY KEY, chrom TEXT, start INTEGER, end INTEGER);
        
        CREATE INDEX idx_snp ON snp(chrom, pos);
        CREATE INDEX idx_res ON restriction_site(chrom, start);
        CREATE INDEX idx_exon ON exon(chrom, start, end);
        CREATE INDEX idx_repeats ON repeats(chrom, start, end);
    """)
    print("✅ 스키마 준비 완료")

def parse_gff3(filename):
    path = os.path.join(RAW_DATA_DIR, filename)
    if not os.path.exists(path): return []
    print(f"📖 Exon 파싱 시작: {filename}")
    data = []
    open_func = gzip.open if filename.endswith('.gz') else open
    try:
        with open_func(path, 'rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split('\t')
                if len(parts) < 9 or parts[2] != 'exon': continue
                chrom, start, end = parts[0], int(parts[3]), int(parts[4])
                attr = parts[8]
                tid = "unknown"
                if "Parent=" in attr:
                    tid = attr.split("Parent=")[1].split(";")[0].replace("transcript:", "")
                data.append((chrom, start, end, tid))
    except Exception as e:
        print(f"❌ Exon 파싱 오류: {e}")
    return data

def parse_vcf(filename):
    path = os.path.join(RAW_DATA_DIR, filename)
    if not os.path.exists(path): return []
    print(f"📖 SNP 파싱 시작: {filename}")
    data = []
    open_func = gzip.open if filename.endswith('.gz') else open
    try:
        with open_func(path, 'rt', encoding='utf-8') as f:
            for line in f:
                if line.startswith("#"): continue
                parts = line.strip().split('\t')
                if len(parts) < 2: continue
                data.append((parts[0], int(parts[1])))
    except Exception as e:
        print(f"❌ SNP 파싱 오류: {e}")
    return data

def parse_repeats_rmsk(filename):
    path = os.path.join(RAW_DATA_DIR, filename)
    if not os.path.exists(path):
        print(f"⚠️  파일 없음: {filename} (Repeats 건너뜀)")
        return []
    
    print(f"📖 Repeats 파싱 시작: {filename}")
    data = []
    open_func = gzip.open if filename.endswith('.gz') else open
    try:
        with open_func(path, 'rt', encoding='utf-8') as f:
            for line in f:
                parts = line.strip().split('\t')
                # UCSC rmsk.txt: col 5(chrom), 6(start), 7(end)
                if len(parts) < 8: continue
                chrom = parts[5]
                start = int(parts[6])
                end = int(parts[7])
                data.append((chrom, start, end))
    except Exception as e:
        print(f"❌ Repeats 파싱 오류: {e}")
    return data

def scan_restriction_sites(fasta_filename):
    path = os.path.join(RAW_DATA_DIR, fasta_filename)
    if not os.path.exists(path):
        print(f"⚠️  파일 없음: {fasta_filename} (제한효소 스캔 건너뜀)")
        return []

    print(f"🕵️ 제한효소 스캔 시작 (FASTA 읽는 중... 시간 소요 예상): {fasta_filename}")
    sites_found = []
    seq_buffer = []
    current_chrom = None
    open_func = gzip.open if fasta_filename.endswith('.gz') else open

    try:
        with open_func(path, 'rt', encoding='utf-8') as f:
            for line in f:
                line = line.strip()
                if not line: continue
                if line.startswith(">"):
                    if current_chrom and seq_buffer:
                        full_seq = "".join(seq_buffer).upper()
                        for name, motif in ENZYMES.items():
                            pos = full_seq.find(motif)
                            while pos != -1:
                                sites_found.append((current_chrom, name, pos+1, pos+len(motif)))
                                pos = full_seq.find(motif, pos+1)
                        print(f"   -> {current_chrom} 스캔 완료")
                        seq_buffer = []
                    current_chrom = line[1:].split()[0]
                else:
                    seq_buffer.append(line)
            
            # 마지막 염색체 처리
            if current_chrom and seq_buffer:
                full_seq = "".join(seq_buffer).upper()
                for name, motif in ENZYMES.items():
                    pos = full_seq.find(motif)
                    while pos != -1:
                        sites_found.append((current_chrom, name, pos+1, pos+len(motif)))
                        pos = full_seq.find(motif, pos+1)
                print(f"   -> {current_chrom} 스캔 완료")

    except Exception as e:
        print(f"❌ 제한효소 스캔 오류: {e}")
    return sites_found

def main():
    conn = get_db_connection()
    init_schema(conn)
    cursor = conn.cursor()

    # 1. Exon (파일명 확인 필)
    exons = parse_gff3("gencode.v49.annotation.gff3.gz") 
    if exons:
        print(f"💾 Exon {len(exons):,}개 저장 중...")
        cursor.executemany("INSERT INTO exon (chrom, start, end, transcript_id) VALUES (?, ?, ?, ?)", exons)
        conn.commit()

    # 2. SNP (파일명 확인 필)
    snps = parse_vcf("clinvar.vcf.gz")
    if snps:
        print(f"💾 SNP {len(snps):,}개 저장 중...")
        cursor.executemany("INSERT INTO snp (chrom, pos) VALUES (?, ?)", snps)
        conn.commit()

    # 3. Repeats (파일명: rmsk.txt.gz)
    repeats = parse_repeats_rmsk("rmsk.txt.gz")
    if repeats:
        print(f"💾 Repeats {len(repeats):,}개 저장 중...")
        cursor.executemany("INSERT INTO repeats (chrom, start, end) VALUES (?, ?, ?)", repeats)
        conn.commit()

    # 4. Restriction Sites (파일명: GRCh38.primary_assembly.genome.fa.gz)
    res_sites = scan_restriction_sites("GRCh38.primary_assembly.genome.fa.gz")
    if res_sites:
        print(f"💾 Restriction Site {len(res_sites):,}개 저장 중...")
        cursor.executemany("INSERT INTO restriction_site (chrom, name, start, end) VALUES (?, ?, ?, ?)", res_sites)
        conn.commit()

    conn.close()
    print(f"\n🎉 최종 DB 구축 완료! 파일 위치: {DB_PATH}")

if __name__ == "__main__":
    main()
