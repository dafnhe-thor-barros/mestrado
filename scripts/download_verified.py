import os
import urllib.request
import pypdf
import io

BIB_DIR = r"C:\ZZ - Antigravity\Mestrado\projeto-mestrado\bibliografia"

def download_and_verify(url, filename, expected_words=None):
    out_path = os.path.join(BIB_DIR, filename)
    if os.path.exists(out_path) and os.path.getsize(out_path) > 10000:
        print(f"[ALREADY EXISTS] {filename} ({os.path.getsize(out_path)} bytes)")
        return True

    print(f"Downloading {filename} from {url}...")
    headers = {
        'User-Agent': 'Mozilla/5.0 (Windows NT 10.0; Win64; x64) AppleWebKit/537.36 (KHTML, like Gecko) Chrome/122.0.0.0 Safari/537.36'
    }
    req = urllib.request.Request(url, headers=headers)
    try:
        with urllib.request.urlopen(req, timeout=30) as resp:
            data = resp.read()
            if len(data) < 5000:
                print(f"[FAIL] Downloaded data too small ({len(data)} bytes)")
                return False
            # Verify PDF
            try:
                reader = pypdf.PdfReader(io.BytesIO(data))
                first_page = reader.pages[0].extract_text()
                print(f"[VERIFIED PDF] {filename} - Pages: {len(reader.pages)} - Size: {len(data)} bytes")
                if expected_words:
                    found = all(w.lower() in first_page.lower() for w in expected_words)
                    if not found:
                        print(f"[WARN] Some expected words not in first page of {filename}")
            except Exception as pe:
                print(f"[FAIL] Not a valid PDF: {pe}")
                return False

            with open(out_path, 'wb') as f:
                f.write(data)
            print(f"[SAVED] {out_path}")
            return True
    except Exception as e:
        print(f"[ERROR] {filename}: {e}")
        return False

def main():
    # 1. Plett 2004
    download_and_verify(
        "http://mocha-java.uccs.edu/dossier/RESEARCH/2004jps1-.pdf",
        "Plett_2004_Extended_Kalman_Filtering_Battery_Management_HEV_LiPB.pdf",
        ["Kalman", "Plett"]
    )

    # 2. Kopp & Lean 2011
    download_and_verify(
        "https://geoffstuff.com/kandl_toa_paper.pdf",
        "Kopp_Lean_2011_New_Lower_Value_Total_Solar_Irradiance_GRL.pdf",
        ["Kopp", "Total Solar Irradiance"]
    )

if __name__ == '__main__':
    main()
