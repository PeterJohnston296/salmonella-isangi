import importlib.util, json, subprocess, sys, tempfile, unittest
from pathlib import Path
ROOT=Path(__file__).resolve().parents[1]
spec=importlib.util.spec_from_file_location("audit",ROOT/"scripts/audit.py")
audit=importlib.util.module_from_spec(spec);spec.loader.exec_module(audit)

class ReleaseTests(unittest.TestCase):
    def test_publication_data(self):
        result=audit.validate()
        self.assertTrue((result.result=="pass").all())
    def test_burden_against_approved_plot(self):
        study,summary,frequency,note=audit.burden()
        counts=dict(zip(summary.group,summary.n))
        self.assertEqual(counts,{"NCBI non-Isangi":810319,"Isangi ST335":224,"Other Isangi":121})
        self.assertEqual(note["wholly_empty_TSV_records_ignored"],12)
        self.assertEqual(note["other_Isangi_unassigned_ST"],6)
    def test_phenix_fasta_extraction(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/"in.fa").write_text(">sampleA\nACGTN\n>reference\nACGTT\n")
            (p/"refs.txt").write_text("reference\n")
            cmd=[sys.executable,str(ROOT/"workflow/lib/fasta_tools.py"),"extract",str(p/"in.fa"),"sampleA",str(p/"refs.txt"),str(p/"sample.fa")]
            subprocess.run(cmd,check=True,capture_output=True,text=True)
            self.assertEqual((p/"sample.fa").read_text(),">sampleA\nACGTN\n")
    def test_fasta_order_validation(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/"in.fa").write_text(">A\nACGTN\n>B\nACGTT\n")
            (p/"ids.txt").write_text("A\nB\n")
            cmd=[sys.executable,str(ROOT/"workflow/lib/fasta_tools.py"),"check",str(p/"in.fa"),str(p/"ids.txt"),str(p/"check.json")]
            subprocess.run(cmd,check=True,capture_output=True,text=True)
            self.assertTrue(json.loads((p/"check.json").read_text())["ordered_membership_verified"])
    def test_invalid_fasta_rejected(self):
        with tempfile.TemporaryDirectory() as tmp:
            p=Path(tmp)
            (p/"in.fa").write_text(">A\nACGTN\n>B\nACGT\n")
            (p/"ids.txt").write_text("A\nB\n")
            cmd=[sys.executable,str(ROOT/"workflow/lib/fasta_tools.py"),"check",str(p/"in.fa"),str(p/"ids.txt"),str(p/"check.json")]
            result=subprocess.run(cmd,capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
    def test_shell_syntax(self):
        for p in (ROOT/"workflow").rglob("*.sh"):
            with self.subTest(script=p.name):
                subprocess.run(["bash","-n",str(p)],check=True,capture_output=True,text=True)
    def test_json_documents(self):
        for p in (ROOT/"data").rglob("*.json"):
            with self.subTest(file=p.name):json.loads(p.read_text())
    def test_no_large_git_files(self):
        for p in ROOT.rglob("*"):
            if p.is_file() and "rebuilt" not in p.parts:
                self.assertLess(p.stat().st_size,100*1024*1024,str(p.relative_to(ROOT)))
    def test_no_sequence_read_or_font_binaries(self):
        for p in ROOT.rglob("*"):
            if p.is_file():
                self.assertNotIn(p.suffix.lower(),[".bam",".cram",".fastq",".fq",".ttf",".otf",".woff",".woff2"],str(p))
if __name__=="__main__":unittest.main()
