"""List / extract members of a remote ZIP via HTTP range requests (no full download).

Usage: python3 -I validation_remote_zip.py <url> list
       python3 -I validation_remote_zip.py <url> extract <member> <out_dir>
Used to pull material/res_test/aggregated.csv (0.8 MB) out of the 192.7 MB AndroWatts zip (Zenodo 14314943).
Caveat: members of a *nested* deflated zip cannot be seeked cheaply; download the outer file instead."""
import sys, io, zipfile, urllib.request, os, ssl

class HttpRangeFile(io.RawIOBase):
    def __init__(self, url):
        self.url = url
        req = urllib.request.Request(url, method="HEAD")
        with urllib.request.urlopen(req, timeout=60) as r:
            self.size = int(r.headers["Content-Length"])
        self.pos = 0
        self.fetched = 0
    def seekable(self): return True
    def readable(self): return True
    def tell(self): return self.pos
    def seek(self, off, whence=0):
        if whence == 0: self.pos = off
        elif whence == 1: self.pos += off
        else: self.pos = self.size + off
        return self.pos
    def read(self, n=-1):
        if n is None or n < 0: n = self.size - self.pos
        if n == 0 or self.pos >= self.size: return b""
        end = min(self.pos + n, self.size) - 1
        req = urllib.request.Request(self.url, headers={"Range": f"bytes={self.pos}-{end}"})
        with urllib.request.urlopen(req, timeout=120) as r:
            data = r.read()
        self.pos += len(data); self.fetched += len(data)
        return data
    def readinto(self, b):
        d = self.read(len(b)); b[:len(d)] = d; return len(d)

def main():
    url = sys.argv[1]; mode = sys.argv[2] if len(sys.argv) > 2 else "list"
    f = HttpRangeFile(url)
    z = zipfile.ZipFile(io.BufferedReader(f, buffer_size=1 << 16))
    if mode == "list":
        for i in z.infolist():
            print(f"{i.file_size/1e6:10.3f} MB  comp {i.compress_size/1e6:10.3f} MB  {i.filename}")
    elif mode == "extract":
        member = sys.argv[3]; outdir = sys.argv[4]
        os.makedirs(outdir, exist_ok=True)
        z.extract(member, outdir)
        print("extracted", member, "to", outdir)
    print(f"[bytes fetched: {f.fetched/1e6:.3f} MB of {f.size/1e6:.1f} MB]", file=sys.stderr)

main()
