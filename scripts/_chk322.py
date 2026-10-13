import gzip
txt = gzip.open(r"C:\Users\hthome\AppData\Local\Temp\bioc322.gz", "rt", errors="replace").read()
for blk in txt.split("\n\n"):
    if blk.startswith("Package: ANCOMBC"):
        for ln in blk.splitlines():
            if ln.startswith(("Package:","Version:","Depends:","Imports:")):
                print(ln.strip()[:200])
