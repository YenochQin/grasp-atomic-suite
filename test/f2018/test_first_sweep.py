"""Run the comparator with controlled executable outputs, including false greens."""
from pathlib import Path
import subprocess
import sys
import tempfile
import unittest

from compare_first_sweep import normalize

PRODUCER = r'''
from pathlib import Path
import struct
name=Path(__file__).name
mode=Path(__file__).parent.name
def record(data): return struct.pack('<i',len(data))+data+struct.pack('<i',len(data))
def wave():
    return record(b'G92RWF')+record(struct.pack('<iidi',1,-1,-1.,1))+record(struct.pack('<3d',0.,1.,0.))+record(struct.pack('<d',0.))
if name=='rwfnestimate':
    Path('rwfn.inp').write_bytes(wave())
elif name=='rangular':
    for n in (30,31,32,33):
        if mode=='missing' and n==33: continue
        sizes=[12,4,8,15,12,12,11,4,8] if n==30 else [15,12,12,15,8]
        if mode=='missing_record': sizes=sizes[:-1]
        Path(f'mcp.{n}').write_bytes(b''.join(record(b'x'*size) for size in sizes))
    Path('rangular.log').write_text('angle 1.0\n')
else:
    Path('rwfn.out').write_bytes(wave())
    Path('rmix.out').write_bytes(record(b'G92MIX')+record(struct.pack('<6i',2,1,1,1,1,1))+
        record(struct.pack('<5i',1,1,1,1,1))+record(struct.pack('<i',1))+
        record(struct.pack('<2d',-1.,0.))+record(struct.pack('<d',1.)))
    Path('rmcdhf.log').write_text('input fixed\n')
    Path('rmcdhf.sum').write_text('Energy '+('2.0' if mode=='energy' else '1.0')+'\n')
    print('Iteration number '+('2' if mode=='iteration' else '1'))
    print('Maximum iterations in SCF Exceeded.')
    if mode=='nan': print('diagnostic NaN')
    if mode=='warning': print('WARNING: changed diagnostic')
    if mode=='truncated': Path('rmix.out').write_bytes(b'bad')
    if mode=='empty': Path('mcp.33').write_bytes(b'')
    if mode=='extra': Path('mcp.34').write_bytes(record(b'fixture'))
print('Wall time:')
print('123 seconds' if mode=='time' else '0 seconds')
print('  Time (Hr/Min/Sec): '+('23/59/59.999' if mode=='time' else '00/00/00.000'))
if mode!='unfinished': print(('RMCDHF' if name.startswith('rmcdhf') else name.upper())+': Execution complete.')
'''


class FirstSweepTests(unittest.TestCase):
    def compare(self, mode: str, *, both: bool = False, optimized: bool = False) -> subprocess.CompletedProcess[str]:
        with tempfile.TemporaryDirectory() as tmp:
            root = Path(tmp)
            dirs = [root / 'reference' / (mode if both else 'good'), root / 'candidate' / mode]
            for directory in dirs:
                directory.mkdir(parents=True)
                for name in ('rwfnestimate', 'rangular', 'rmcdhf', 'rmcdhf_mem'):
                    binary = directory / name
                    binary.write_text('#!'+sys.executable+'\n'+PRODUCER)
                    binary.chmod(0o700)
            result = subprocess.run([sys.executable, *(['-O'] if optimized else []),
                str(Path(__file__).with_name('compare_first_sweep.py')), '--reference-bindir', str(dirs[0]),
                '--candidate-bindir', str(dirs[1]), '--output', str(root/'result')], capture_output=True, text=True)
            self.assertTrue((root/'result/comparison.json').is_file())
            return result

    def test_only_documented_metadata_can_differ(self) -> None:
        for mode in ('good', 'time'):
            result = self.compare(mode)
            self.assertEqual(result.returncode, 0, result.stderr)
        self.assertNotEqual(normalize('time step 1\n'), normalize('time step 2\n'))

    def test_changed_science_warning_iteration_or_completion_fails(self) -> None:
        for mode in ('warning', 'energy', 'iteration', 'unfinished'):
            with self.subTest(mode=mode):
                self.assertNotEqual(self.compare(mode).returncode, 0)

    def test_matching_invalid_or_incomplete_outputs_fail(self) -> None:
        for mode in ('nan', 'missing', 'empty', 'extra', 'truncated', 'unfinished', 'missing_record'):
            with self.subTest(mode=mode):
                self.assertNotEqual(self.compare(mode, both=True, optimized=True).returncode, 0)


if __name__ == '__main__':
    unittest.main()
