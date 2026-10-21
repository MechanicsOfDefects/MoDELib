"""Build in an isolated directory so historical auxiliary files cannot interfere."""
from pathlib import Path
import argparse, shutil, subprocess, tempfile
p=Path(__file__).resolve().parents[1]
a=argparse.ArgumentParser();a.add_argument('--publish',action='store_true',help='copy successful PDF to manual_beta');args=a.parse_args()
b=Path(tempfile.mkdtemp(prefix='modelib-manual-'))
shutil.copy2(p/'MoDELib_manual.tex',b)
shutil.copy2(p/'manual_cited.bib',b)
for d in ['tex','figures']:shutil.copytree(p/d,b/d)
log=b/'build-output.txt'
with log.open('w') as f:
    result=subprocess.run(['latexmk','-pdf','-interaction=nonstopmode','-halt-on-error','MoDELib_manual.tex'],cwd=b,stdout=f,stderr=subprocess.STDOUT)
print('Build directory:',b)
print('Build log:',log)
if result.returncode:raise SystemExit(result.returncode)
if args.publish:shutil.copy2(b/'MoDELib_manual.pdf',p/'MoDELib_manual.pdf')
print('PDF:',p/'MoDELib_manual.pdf' if args.publish else b/'MoDELib_manual.pdf')
