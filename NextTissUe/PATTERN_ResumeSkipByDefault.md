# Pattern: Resume/Continue with Skip-By-Default (Safe Resumable Batch Jobs)

**Last Updated:** 2026-01-01
**Source:** NextTissUe/run_parameter_sweep.py
**Maturity:** Production-ready, tested on 32-job parameter sweeps with interruptions

---

## Overview

This pattern implements a **safety-first resumable batch job system** that automatically skips completed work by default, preventing accidental re-execution of expensive computations. It uses simple file markers for completion tracking and maintains an append-only audit log.

### Core Philosophy: **Safe By Default**

Traditional approach (dangerous):
```bash
# Default: Re-runs everything
./batch_job.py --all

# Need flag to skip completed work
./batch_job.py --all --resume
```

**This pattern (safe):**
```bash
# Default: Automatically skips completed work
./batch_job.py --all

# Need explicit flag to force re-run
./batch_job.py --all --force-rerun
```

**Why invert the logic?**
- ✅ Prevents accidental data loss from expensive simulations
- ✅ Makes resume the common case (just re-run same command)
- ✅ Requires explicit intent to destroy completed work
- ✅ Matches user expectations for long-running batch jobs

### Key Features

- ✅ **Atomic completion markers** (simple file existence check)
- ✅ **Skip-by-default behavior** (safety first)
- ✅ **Explicit force-rerun override** (conscious destruction)
- ✅ **Append-only audit log** (CSV manifest)
- ✅ **Interrupt-safe** (each job is independently resumable)
- ✅ **No partial state** (completion marker only after all artifacts saved)

### When to Use

✅ **Perfect for:**
- Scientific simulations (hours per job)
- ML training runs (expensive to repeat)
- Video rendering, data processing
- Multi-stage pipelines with checkpoints
- Any batch job where interruption is likely

❌ **Not suitable for:**
- Jobs completing in seconds (overhead not worth it)
- Truly stateless operations (no output to preserve)
- Strictly ordered dependencies (would need DAG tracking)
- Concurrent execution without locking

---

## Complete Implementation

### Minimal Reusable Pattern

```python
from pathlib import Path
from datetime import datetime
import csv

class ResumableBatchRunner:
    """Universal pattern for skip-by-default resumable batch processing

    This class provides the core resume functionality that can be adapted
    to any batch job system.
    """

    def __init__(self, output_base, manifest_name='batch_manifest.csv'):
        """Initialize batch runner

        Args:
            output_base: Base directory for job outputs (Path or str)
            manifest_name: CSV manifest filename
        """
        self.output_base = Path(output_base)
        self.manifest_file = self.output_base / manifest_name

    def is_complete(self, job_dir):
        """Check if job has completion marker

        Args:
            job_dir: Job directory path (Path or str)

        Returns:
            True if COMPLETE marker exists, False otherwise

        Why a file instead of database/state file?
        - Simple: No dependencies, works everywhere
        - Atomic: File existence is atomic on most filesystems
        - Visible: User can see what's done with 'ls'
        - Robust: Survives process crashes, network failures
        """
        return (Path(job_dir) / 'COMPLETE').exists()

    def mark_complete(self, job_dir):
        """Create completion marker file

        Args:
            job_dir: Job directory path (Path or str)

        CRITICAL: Only call this AFTER all job artifacts are safely saved.
        The marker should be the last thing written.
        """
        marker = Path(job_dir) / 'COMPLETE'
        with open(marker, 'w') as f:
            f.write(f"Completed at {datetime.now().isoformat()}\n")

    def update_manifest(self, job_params, status, start_time=None, end_time=None, duration=None, job_dir=None):
        """Append job record to CSV manifest

        Args:
            job_params: Dict of job parameters
            status: 'completed', 'failed', 'error'
            start_time: datetime or None
            end_time: datetime or None
            duration: float (seconds) or None
            job_dir: Job directory path

        Why append-only?
        - Survives crashes (partial writes still visible)
        - No need for locking (safe for sequential writes)
        - Provides audit trail (see when job was attempted)
        - Simple: No database, just CSV

        Why no deduplication?
        - Skipped jobs are re-logged as 'completed' (intentional)
        - Shows when sweep was resumed in audit trail
        - CSV tools can filter duplicates if needed
        """
        file_exists = self.manifest_file.exists()

        row = {
            **job_params,  # Include all job parameters
            'status': status,
            'start_time': start_time.isoformat() if start_time else '',
            'end_time': end_time.isoformat() if end_time else '',
            'duration_s': f"{duration:.1f}" if duration else '',
            'output_dir': str(job_dir) if job_dir else ''
        }

        # Ensure output directory exists
        self.output_base.mkdir(parents=True, exist_ok=True)

        # Append to CSV (create with header if new file)
        with open(self.manifest_file, 'a', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=row.keys())
            if not file_exists:
                writer.writeheader()
            writer.writerow(row)

    def execute_batch(self, jobs, run_job_fn, force_rerun=False):
        """Execute batch of jobs with skip-by-default logic

        Args:
            jobs: Iterable of job dictionaries
            run_job_fn: Function(job_dict) -> (success: bool, job_dir: Path)
            force_rerun: If True, ignore COMPLETE markers and re-run everything

        Returns:
            Dict with statistics: {'completed': N, 'failed': M, 'skipped': K}

        Core Logic:
            1. Check COMPLETE marker
            2. If exists and not force_rerun: skip and log
            3. Otherwise: run job, mark complete if successful
            4. Append to manifest regardless
        """
        completed = 0
        failed = 0
        skipped = 0

        for job in jobs:
            job_dir = self._get_job_dir(job)

            # CORE PATTERN: Skip-by-default logic
            if not force_rerun and self.is_complete(job_dir):
                print(f"SKIPPED: {job_dir.name} (use --force-rerun to override)")
                self.update_manifest(job, 'completed', job_dir=job_dir)
                skipped += 1
                completed += 1  # Count toward overall progress
                continue

            # Execute job
            start_time = datetime.now()
            try:
                success, actual_job_dir = run_job_fn(job)
                end_time = datetime.now()
                duration = (end_time - start_time).total_seconds()

                if success:
                    self.mark_complete(actual_job_dir)
                    self.update_manifest(job, 'completed', start_time, end_time, duration, actual_job_dir)
                    completed += 1
                else:
                    self.update_manifest(job, 'failed', start_time, end_time, duration, actual_job_dir)
                    failed += 1

            except Exception as e:
                end_time = datetime.now()
                duration = (end_time - start_time).total_seconds()
                self.update_manifest(job, 'error', start_time, end_time, duration, job_dir)
                print(f"ERROR: {job_dir.name}: {e}")
                failed += 1

        return {'completed': completed, 'failed': failed, 'skipped': skipped}

    def _get_job_dir(self, job):
        """Override this to customize job directory naming"""
        # Example: use job ID or hash of parameters
        job_id = job.get('id', hash(frozenset(job.items())))
        return self.output_base / f"job_{job_id}"
```

---

## Design Evolution: Why Skip-by-Default?

### Iteration 1: No Resume Capability (Original Problem)

```python
# PROBLEM: No way to resume interrupted batch
for job in jobs:
    run_job(job)  # Crash at job 15/32 = lose all progress
```

**Issues:**
- ❌ Any interruption (crash, Ctrl+C, power loss) loses all work
- ❌ Must manually track which jobs completed
- ❌ Easy to accidentally re-run expensive jobs

### Iteration 2: Opt-In Resume (Traditional Approach)

```python
# DEFAULT: Re-run everything
./batch.py --all

# RESUME: Need to remember flag every time
./batch.py --all --resume
```

**Code:**
```python
if args.resume and is_complete(job_dir):
    skip_job()
else:
    run_job()
```

**Issues:**
- ⚠️ Easy to forget `--resume` flag
- ⚠️ Accidental re-runs waste hours/days of compute
- ⚠️ Resume should be the common case, not special case

### Iteration 3: Skip-by-Default (Final Solution)

```python
# DEFAULT: Automatically skip completed
./batch.py --all

# FORCE RERUN: Explicit intent required
./batch.py --all --force-rerun
```

**Code:**
```python
if not force_rerun and is_complete(job_dir):
    skip_job()  # DEFAULT: Skip
else:
    run_job()   # Only if forced or incomplete
```

**Why this is better:**
- ✅ **Safe by default** - never waste compute accidentally
- ✅ **Resume is implicit** - just re-run same command
- ✅ **Explicit destruction** - `--force-rerun` makes intent clear
- ✅ **Matches expectations** - users expect long jobs to resume

**The Naming Inversion:**
```python
# Parameter is named 'force_rerun' for clarity
def execute_batch(self, jobs, force_rerun=False):
    if not force_rerun and self.is_complete(job_dir):
        # When force_rerun=False (default): SKIP completed
        # When force_rerun=True (explicit): RE-RUN everything
```

Old code used `resume` parameter but with inverted logic - confusing. New pattern uses `force_rerun` for clarity.

---

## Technical Deep Dive

### COMPLETE Marker System

**Why a file instead of database?**

| Approach | Pros | Cons |
|----------|------|------|
| **File marker** | Simple, visible, atomic, no dependencies | Not queryable |
| Database | Queryable, centralized | Requires setup, single point of failure |
| State file | Structured data | Requires parsing, locking, can corrupt |

**File marker wins for:**
- Simplicity (works everywhere)
- Atomicity (file existence is atomic)
- Visibility (`ls` shows what's done)
- Robustness (survives crashes)

**Marker Content:**
```
Completed at 2026-01-01T14:32:15.123456
```

**Could be enhanced with:**
```
Completed at 2026-01-01T14:32:15.123456
Duration: 1234.5 seconds
Version: 1.2.3
Checksum: a1b2c3d4...
```

But minimal content is sufficient for most use cases.

**Atomicity Guarantee:**

File creation is atomic on most filesystems:
```python
with open(marker, 'w') as f:
    f.write(timestamp)
# File appears atomically after close()
```

**Critical: Marker Last**
```python
# CORRECT order:
run_simulation()
save_all_outputs()
verify_outputs_exist()
mark_complete()  # LAST STEP

# WRONG order:
mark_complete()  # Marked complete but no data!
run_simulation()
```

### Manifest Tracking Strategy

**CSV Format:**
```csv
job_id,param1,param2,status,start_time,end_time,duration_s,output_dir
job_1,500,10,completed,2026-01-01T10:00:00,2026-01-01T10:03:00,180.0,/path/to/job_1
job_2,500,20,failed,2026-01-01T10:03:00,2026-01-01T10:05:30,150.0,/path/to/job_2
job_1,500,10,completed,,,/path/to/job_1
```

**Note**: `job_1` appears twice - once from original run, once from resume (skipped).

**Why append-only?**
- No need to update existing rows
- Audit trail of all attempts
- Survives crashes mid-write
- Simple implementation (no locking)

**Deduplication in analysis:**
```python
import pandas as pd
df = pd.read_csv('manifest.csv')
# Get latest status for each job
latest = df.groupby('job_id').tail(1)
```

### Progress Counter Integration

**The Dual Counter Pattern:**

```python
completed = 0  # Jobs that finished successfully
skipped = 0    # Jobs that were skipped (already complete)
failed = 0     # Jobs that failed/errored

# Skipped jobs increment BOTH counters
if not force_rerun and is_complete(job_dir):
    skipped += 1
    completed += 1  # Also count as completed for progress
    continue

# Overall progress calculation
runs_done = completed + failed  # All jobs attempted
total_progress = runs_done / total_jobs
```

**Why increment both `skipped` and `completed`?**
- Skipped jobs are "done" from the user's perspective
- They count toward overall batch progress
- But track separately for reporting

**Output example:**
```
Overall: 15/32 (C:12 F:1 S:2)
         ^^^^  ^^^ ^^ ^^^
         done  ok  bad already done
```

### ETA Calculation Interaction

**Skipped jobs affect ETA calculation:**

```python
# Collect durations from EXECUTED jobs only
run_durations = []  # Does NOT include skipped jobs

for job in jobs:
    if skipped:
        # Don't add to run_durations
        continue

    start = time.time()
    run_job()
    run_durations.append(time.time() - start)

# ETA based on actual execution times
if run_durations:
    avg_duration = sum(run_durations) / len(run_durations)
    remaining_jobs = total_jobs - len(run_durations) - skipped
    eta_seconds = remaining_jobs * avg_duration
```

**Critical:** Skipped jobs have zero duration but count toward progress.

**Edge case: All skipped**
```python
if not run_durations:
    # No timing data yet
    eta = "calculating..."
else:
    eta = estimate_from(run_durations)
```

---

## Edge Cases and Solutions

### 1. Crash During Job Execution

**Scenario:** Process killed while job is running

**Behavior:**
```python
start_job()
# Process killed here - no COMPLETE marker created
mark_complete()  # Never reached
```

**Result:** Job will be retried on resume (correct behavior)

**Why this works:** Marker is written LAST, after all outputs. No marker = incomplete = retry.

### 2. Crash After Job Completes

**Scenario:** Process killed after job finishes but before manifest update

**Behavior:**
```python
run_job()
mark_complete()  # Marker created
# Process killed here
update_manifest()  # Never reached
```

**Result:** Job marked complete, skipped on resume, manifest updated on skip

**Why this works:** Marker is source of truth. Manifest is audit log, not state.

### 3. Interrupted Resume

**Scenario:** Resume operation itself is interrupted

**Behavior:**
- Jobs completed before interruption: have markers, will be skipped
- Jobs interrupted during execution: no markers, will retry
- Jobs not yet started: will run in next resume

**Why this works:** Each job is independently tracked. No global state.

### 4. Concurrent Execution (Not Supported)

**Scenario:** Two processes run same batch simultaneously

**Problem:**
```python
# Process A
if not is_complete(job_dir):  # Check
    # Process B checks here too - both see incomplete
    run_job()  # Both run same job!
    mark_complete()
```

**Result:** Both processes run same job, last one overwrites

**Solution:** This pattern is NOT designed for concurrent execution. Use:
- Distributed lock (Redis, file lock)
- Job queue (Celery, RQ)
- Unique job dirs per process

### 5. Filesystem Delays

**Scenario:** Network filesystem delayed sync

**Problem:**
```python
# Process on machine A
mark_complete(job_dir)

# Process on machine B (microseconds later)
if is_complete(job_dir):  # File not visible yet
    run_job()  # Duplicate run!
```

**Mitigation:**
```python
import time

def mark_complete(job_dir):
    marker = job_dir / 'COMPLETE'
    with open(marker, 'w') as f:
        f.write(...)
    # Force sync
    f.flush()
    os.fsync(f.fileno())

    # Verify marker is readable
    time.sleep(0.1)
    assert marker.exists()
```

### 6. Partial Output Files

**Scenario:** Job writes partial data but exits with code 0

**Problem:**
```python
run_job()  # Writes partial output, exits 0
if returncode == 0:
    mark_complete()  # Marked complete but data incomplete!
```

**Solution:** Validate outputs before marking complete
```python
success = run_job()
if success and validate_outputs(job_dir):
    mark_complete()
else:
    # Don't mark complete - will retry
    log_error("Validation failed")
```

---

## Universal Pattern for Any Batch System

This pattern abstracts to handle ANY long-running batch job:

```python
class GenericBatchSystem:
    """Universal resumable batch pattern - adapt to your needs"""

    def __init__(self, output_base):
        self.output_base = Path(output_base)
        self.manifest = self.output_base / 'manifest.csv'

    def run_batch(self, jobs, execute_fn, force_rerun=False):
        """Main execution loop

        Args:
            jobs: List of job specifications (dicts, objects, etc.)
            execute_fn: Callable that runs one job
            force_rerun: Skip safety check
        """
        for job in jobs:
            job_dir = self._prepare_job_dir(job)

            # CORE PATTERN: Skip-by-default
            if not force_rerun and self._is_complete(job_dir):
                self._log_skip(job)
                continue

            # Execute
            try:
                execute_fn(job, job_dir)
                self._mark_complete(job_dir)
                self._log_success(job, job_dir)
            except Exception as e:
                self._log_failure(job, job_dir, e)

    def _is_complete(self, job_dir):
        return (job_dir / 'COMPLETE').exists()

    def _mark_complete(self, job_dir):
        (job_dir / 'COMPLETE').write_text(f"{datetime.now().isoformat()}\n")

    def _prepare_job_dir(self, job):
        """Create and return job directory - customize naming here"""
        job_dir = self.output_base / self._get_job_name(job)
        job_dir.mkdir(parents=True, exist_ok=True)
        return job_dir

    def _get_job_name(self, job):
        """Override this for custom naming logic"""
        raise NotImplementedError

    def _log_skip(self, job):
        """Append skip record to manifest"""
        pass

    def _log_success(self, job, job_dir):
        """Append success record to manifest"""
        pass

    def _log_failure(self, job, job_dir, error):
        """Append failure record to manifest"""
        pass
```

**Adapt by implementing:**
1. `_get_job_name()` - How to name job directories
2. `_log_*()` methods - How to track progress
3. `execute_fn` - What to actually run

---

## Usage Examples

### Example 1: Parameter Sweep

```python
from batch_runner import ResumableBatchRunner

# Define parameter grid
jobs = [
    {'param_a': a, 'param_b': b, 'param_c': c}
    for a in [1, 2, 3]
    for b in [10, 20, 30]
    for c in [100, 200]
]

def run_simulation(job):
    """Run one simulation with given parameters"""
    job_dir = Path(f"run_a{job['param_a']}_b{job['param_b']}_c{job['param_c']}")
    job_dir.mkdir(exist_ok=True)

    # Run simulation
    cmd = ['./simulate', '-a', str(job['param_a']), '-b', str(job['param_b'])]
    result = subprocess.run(cmd, cwd=job_dir)

    return (result.returncode == 0, job_dir)

# Execute batch
runner = ResumableBatchRunner('output/')
stats = runner.execute_batch(jobs, run_simulation, force_rerun=False)
print(f"Completed: {stats['completed']}, Failed: {stats['failed']}, Skipped: {stats['skipped']}")
```

### Example 2: ML Training Runs

```python
experiments = [
    {'lr': 0.001, 'batch_size': 32, 'epochs': 100},
    {'lr': 0.01, 'batch_size': 64, 'epochs': 100},
    # ... more configs
]

def train_model(experiment):
    exp_dir = Path(f"exp_lr{experiment['lr']}_bs{experiment['batch_size']}")
    exp_dir.mkdir(exist_ok=True)

    # Train
    model = create_model()
    train(model, lr=experiment['lr'], batch_size=experiment['batch_size'])

    # Save
    model.save(exp_dir / 'model.pth')
    save_metrics(exp_dir / 'metrics.json')

    return (True, exp_dir)

runner = ResumableBatchRunner('experiments/')
runner.execute_batch(experiments, train_model)
```

### Example 3: Data Processing Pipeline

```python
files_to_process = ['data1.csv', 'data2.csv', 'data3.csv']
jobs = [{'input': f, 'output': f.replace('.csv', '_processed.csv')} for f in files_to_process]

def process_file(job):
    job_dir = Path('processed') / Path(job['input']).stem
    job_dir.mkdir(parents=True, exist_ok=True)

    # Process
    df = pd.read_csv(job['input'])
    result = expensive_processing(df)
    result.to_csv(job_dir / 'output.csv')

    return (True, job_dir)

runner = ResumableBatchRunner('processed/')
runner.execute_batch(jobs, process_file)
```

---

## Integration with Progress Display

Combine with two-line progress bar pattern:

```python
from progress_display import LiveProgressDisplay
from batch_runner import ResumableBatchRunner

display = LiveProgressDisplay()
runner = ResumableBatchRunner('output/')

for i, job in enumerate(jobs):
    # Check if skippable
    job_dir = runner._get_job_dir(job)
    if not force_rerun and runner.is_complete(job_dir):
        # Update progress but don't run
        line1 = f"Job {i+1}/{len(jobs)}: {job_dir.name}"
        line2 = f"SKIPPED (already complete)"
        display.update(line1, line2)
        continue

    # Run with progress monitoring
    line1 = f"Job {i+1}/{len(jobs)}: Running {job_dir.name}..."
    display.update(line1, f"Status: In progress")

    success, actual_dir = run_job(job)

    if success:
        runner.mark_complete(actual_dir)

display.finish()
```

---

## Common Pitfalls

### ❌ Marking Complete Too Early

```python
# BAD: Marked complete before outputs saved
def run_job(job):
    mark_complete(job_dir)  # WRONG - nothing saved yet!
    save_outputs()

# GOOD: Marker is last step
def run_job(job):
    save_outputs()
    verify_outputs()
    mark_complete(job_dir)  # LAST
```

### ❌ Not Handling Exceptions

```python
# BAD: Exception prevents cleanup, partial state
def run_job(job):
    start_processing()
    # Exception here = no marker, but partial files exist!
    mark_complete()

# GOOD: Clean up on failure
def run_job(job):
    try:
        start_processing()
        mark_complete()
    except Exception:
        cleanup_partial_files()
        raise
```

### ❌ Modifying Manifest Instead of Appending

```python
# BAD: Re-writing entire file (can corrupt)
def update_manifest(job, status):
    data = read_csv('manifest.csv')
    data.append({'job': job, 'status': status})
    write_csv('manifest.csv', data)  # Can corrupt on crash

# GOOD: Append-only (crash-safe)
def update_manifest(job, status):
    with open('manifest.csv', 'a') as f:
        writer = csv.writer(f)
        writer.writerow([job, status])
```

### ❌ Trusting Exit Codes Blindly

```python
# BAD: Exit code 0 doesn't guarantee success
result = subprocess.run(cmd)
if result.returncode == 0:
    mark_complete()  # But outputs might be corrupt!

# GOOD: Validate outputs
result = subprocess.run(cmd)
if result.returncode == 0 and validate_outputs(job_dir):
    mark_complete()
```

---

## Command-Line Interface Pattern

Standard CLI for resumable batch scripts:

```python
import argparse

parser = argparse.ArgumentParser()
parser.add_argument('--all', action='store_true',
                   help='Run all jobs (skips completed by default)')
parser.add_argument('--force-rerun', action='store_true',
                   help='Force re-run even if COMPLETE marker exists')
parser.add_argument('--job-id', type=str,
                   help='Run specific job only')

args = parser.parse_args()

runner = ResumableBatchRunner('output/')

if args.job_id:
    # Run single job
    job = get_job_by_id(args.job_id)
    runner.execute_batch([job], run_fn, force_rerun=args.force_rerun)
elif args.all:
    # Run all jobs
    runner.execute_batch(all_jobs, run_fn, force_rerun=args.force_rerun)
```

**Usage:**
```bash
# Normal: skip completed
./batch.py --all

# Force rerun all
./batch.py --all --force-rerun

# Rerun one specific job
./batch.py --job-id exp_123 --force-rerun
```

---

## Summary

This pattern provides **production-ready resumable batch execution** with:

- ✅ **Skip-by-default** - Never accidentally waste compute
- ✅ **Simple file markers** - Atomic completion tracking
- ✅ **Append-only manifest** - Crash-safe audit log
- ✅ **Interrupt-safe** - Resume from any point
- ✅ **Explicit override** - Require flag to force re-run
- ✅ **Independent jobs** - No global state

**Key Design Principle:** Make the safe behavior the default, require explicit action for dangerous operations.

**Copy this documentation to any project** with long-running batch jobs to add robust resume capability.

**Critical Insight:** The logic inversion from opt-in resume to skip-by-default is the core innovation. It aligns system behavior with user expectations and prevents costly mistakes.
