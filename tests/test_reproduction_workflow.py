"""Check the real shell driver's dependency graph without submitting jobs."""
import json
import os
from pathlib import Path
import shlex
import subprocess
import tempfile
import unittest

ROOT = Path(__file__).resolve().parents[1]

class WorkflowTests(unittest.TestCase):
    def test_only_unknown_solver_status_is_retried(self):
        with tempfile.TemporaryDirectory(dir=os.environ.get('TMPDIR', str(Path.home()/'tmp'))) as directory:
            base=Path(directory)
            runner=base/'fake-python'
            runner.write_text('''#!/usr/bin/env python3
import os, pathlib, sys
budget=sys.argv[sys.argv.index('--solver-seconds')+1]
p=pathlib.Path(os.environ['ATTEMPTS'])
with p.open('a') as f: f.write(budget+'\\n')
if budget=='120':
    print(os.environ['FAILURE'])
    sys.exit(1)
''')
            runner.chmod(0o755)
            tasks=base/'tasks'; tasks.write_text('deepgo-se cc 7955\n')
            config=base/'config'
            values={'CODE_DIR':str(ROOT),'DATA_DIR':str(base),'RESOURCE_DIR':str(base),
                    'OUT_DIR':str(base/'results'),'PYTHON':str(runner),'TASKS':str(tasks)}
            config.write_text(''.join(k+'='+shlex.quote(v)+'\n' for k,v in values.items()))
            attempts=base/'attempts'
            env={**os.environ,'AUDIT_CONFIG':str(config),'SLURM_JOB_ID':'test','SLURM_ARRAY_TASK_ID':'0','ATTEMPTS':str(attempts)}
            command=['bash',str(ROOT/'slurm/revision_audit.slurm')]
            unknown='RuntimeError: CP-SAT did not produce a valid repair: UNKNOWN'
            result=subprocess.run(command,env={**env,'FAILURE':unknown},capture_output=True,text=True)
            self.assertEqual(result.returncode,0,result.stderr)
            self.assertEqual(attempts.read_text().splitlines(),['120','900'])
            attempts.unlink()
            result=subprocess.run(command,env={**env,'FAILURE':'ValueError: invalid inputs'},capture_output=True,text=True)
            self.assertNotEqual(result.returncode,0)
            self.assertEqual(attempts.read_text().splitlines(),['120'])

    def test_generation_depends_on_every_required_analysis(self):
        with tempfile.TemporaryDirectory(prefix='workflow test ', dir=os.environ.get('TMPDIR', str(Path.home()/'tmp'))) as directory:
            base = Path(directory)
            executable = base/'sbatch'
            executable.write_text('''#!/usr/bin/env python3
import json, os, pathlib, sys
p=pathlib.Path(os.environ['MOCK_SUBMISSIONS'])
rows=json.loads(p.read_text()) if p.exists() else []
rows.append(sys.argv[1:])
p.write_text(json.dumps(rows))
print(1000+len(rows))
''')
            executable.chmod(0o755)
            fixture=base/'required.json'; fixture.write_text('{}')
            (base/'reference-manifest.json').write_text('{}')
            config=base/'config.env'
            values={'CODE_DIR':str(ROOT),'DATA_DIR':str(base),'RESOURCE_DIR':str(base),
                'OUT_DIR':str(base/'results'),'PYTHON':'python3','GAEF_DIR':str(base),
                'FUNTAXIS_MANIFEST':str(fixture),'PROTEOME_ROOT':str(base),'PROTEOME_INDEX':str(fixture),
                'CP_SUBUNIT_INDEX':str(fixture),'ANNOTATION_DIR':str(base),
                'SLURM_ACCOUNT':'test','SLURM_PARTITION':'test','AUDIT_WORKERS':'2','AUDIT_MEMORY':'48G'}
            config.write_text(''.join(k+'='+shlex.quote(v)+'\n' for k,v in values.items()))
            log=base/'submissions.json'
            env={**os.environ,'PATH':str(base)+os.pathsep+os.environ['PATH'],'MOCK_SUBMISSIONS':str(log)}
            subprocess.run(['bash',str(ROOT/'reproducibility/run_all_slurm.sh'),str(config)],env=env,check=True,capture_output=True,text=True)
            calls=json.loads(log.read_text())
            self.assertEqual([call[-1] for call in calls],['workers','evidence','aggregate','lineage','membership','validate-exports','generate'])
            self.assertIn('--mem=48G',calls[0])
            self.assertIn('--dependency=afterok:1001',calls[2])
            self.assertIn('--dependency=afterok:1001',calls[3])
            self.assertIn('--dependency=afterok:1001:1002',calls[4])
            self.assertIn('--dependency=afterok:1001',calls[5])
            self.assertIn('--dependency=afterok:1003:1004:1005:1002:1006',calls[6])
            self.assertEqual(len((base/'results/jobs.txt').read_text().splitlines()),7)
            rerun=subprocess.run(['bash',str(ROOT/'reproducibility/run_all_slurm.sh'),str(config)],env=env,capture_output=True,text=True)
            self.assertNotEqual(rerun.returncode,0)
            self.assertEqual(len(json.loads(log.read_text())),7)

if __name__ == '__main__':
    unittest.main()
