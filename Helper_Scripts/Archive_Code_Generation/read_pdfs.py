import PyPDF2
import re
import sys

def extract_pdf(filename, patterns):
    try:
        with open(filename, 'rb') as f:
            reader = PyPDF2.PdfReader(f)
            text = "\n".join([p.extract_text() for p in reader.pages if p.extract_text()])
            for pattern in patterns:
                print(f"\nSearching for '{pattern}' in {filename}:")
                matches = list(re.finditer(pattern, text, re.IGNORECASE | re.DOTALL))
                if not matches:
                    print("No matches found.")
                for m in matches[:5]: # limit to first 5 matches per pattern
                    start = max(0, m.start() - 300)
                    end = min(len(text), m.end() + 300)
                    print(f"--- MATCH ---")
                    print(text[start:end])
    except Exception as e:
        print(f"Error reading {filename}: {e}")

extract_pdf("spin_squeezing.pdf", [r"Holland", r"HB state", r"Burnett"])
extract_pdf("Humphreys et al. - 2013 - Quantum Enhanced Multiple Phase Estimation.pdf", [r"equation \(10", r"equation 10", r"alpha", r"basis\("])
