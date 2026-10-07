"""Draw every non-sterol lipid in a collection, and report what does not work."""
import glob, os, sys, traceback
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import lipid2d

root, out = sys.argv[1], sys.argv[2]
os.makedirs(out, exist_ok=True)
ok, fails = [], {}
for d in sorted(glob.glob(os.path.join(root, '*'))):
    if not os.path.isdir(d):
        continue
    e = os.path.basename(d)
    pdbs = sorted(p for p in glob.glob(os.path.join(d, '*.pdb')) if not p.endswith('-init.pdb'))
    if not pdbs:
        continue
    try:
        n, c, h = lipid2d.draw(pdbs[0], os.path.join(out, f'{e}.svg'))
        ok.append((e, n, c, h))
    except Exception as exc:
        fails.setdefault(f'{type(exc).__name__}: {str(exc).split(chr(10))[0][:70]}', []).append(e)
print(f'drew {len(ok)}, failed {sum(len(v) for v in fails.values())}')
for msg, entries in sorted(fails.items(), key=lambda kv: -len(kv[1])):
    print(f'  [{len(entries):3d}] {msg}')
    print(f'        {", ".join(sorted(entries)[:14])}{" ..." if len(entries) > 14 else ""}')
import collections
print('chains drawn:', dict(sorted(collections.Counter(c for _, _, c, _ in ok).items())))
