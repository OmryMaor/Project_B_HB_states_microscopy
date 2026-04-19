import PyPDF2
import os

pdf_dir = r"c:\Users\omrym\Technion\Semester9\Project B\HB_states_ProjectB\Project_A_files\Report_And_Presentation"
out_dir = r"c:\Users\omrym\Technion\Semester9\Project B\HB_states_ProjectB\Helper_Scripts\Archive_Code_Generation"

pdfs = ["Final Presentation.pdf", "Final Report.pdf", "הכנת מצגת ודוח פרויקט.pdf"]

for pdf in pdfs:
    path = os.path.join(pdf_dir, pdf)
    out_path = os.path.join(out_dir, pdf.replace('.pdf', '.txt'))
    try:
        with open(path, 'rb') as f:
            reader = PyPDF2.PdfReader(f)
            text = []
            for i, page in enumerate(reader.pages):
                text.append(f"--- PAGE {i+1} ---")
                text.append(page.extract_text() or "")
            with open(out_path, 'w', encoding='utf-8') as out_f:
                out_f.write("\n".join(text))
        print(f"Dumped {pdf}")
    except Exception as e:
        print(f"Failed {pdf}: {e}")
