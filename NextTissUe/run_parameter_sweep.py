#!/usr/bin/env python3
"""
run_parameter_sweep.py - Automate multi-parameter LAMMPS sweeps for NextTissUe

This script runs sequential LAMMPS simulations with different parameter combinations
(karea, pre_fac, fd) organized in structured directories with resume capability.
"""

import argparse
import subprocess
import os
import sys
from pathlib import Path
import csv
from datetime import datetime, timedelta
import shutil
import re
import itertools
import yaml

class ParameterSweep:
    def __init__(self, campaigns, lammps_exe, np, input_file, output_base,
                 main_run_steps, dump_interval):
        """Initialize sweep configuration with campaign definitions"""
        self.campaigns = campaigns
        self.lammps_exe = lammps_exe
        self.np = np
        self.input_file = input_file
        self.output_base = Path(output_base)
        self.main_run_steps = main_run_steps
        self.dump_interval = dump_interval
        self.nexttissue_dir = Path(__file__).parent
        self.manifest_file = self.output_base / "sweep_manifest.csv"
        self.log_file = self.nexttissue_dir / "sweep.log"

    def generate_run_combinations(self, campaign_id):
        """Generate all parameter combinations for a campaign

        Returns: List of dicts like:
          [{'karea': 500, 'pre_fac': 10, 'fd': 200},
           {'karea': 500, 'pre_fac': 10, 'fd': 500}, ...]
        """
        campaign = self.campaigns[campaign_id]
        fixed_params = campaign['fixed']
        sweep_params = campaign['sweep']

        # Get all parameter names and values
        all_params = {}
        all_params.update(fixed_params)

        # Generate cartesian product of sweep parameters
        sweep_keys = list(sweep_params.keys())
        sweep_values = [sweep_params[k] for k in sweep_keys]

        combinations = []
        for values in itertools.product(*sweep_values):
            params = dict(all_params)  # Start with fixed params
            for key, val in zip(sweep_keys, values):
                params[key] = val
            params['campaign'] = campaign_id
            combinations.append(params)

        return combinations

    def get_run_dir_name(self, params):
        """Generate directory name: karea{val}_pre_fac{val}_fd{val}"""
        return f"karea{int(params['karea'])}_pre_fac{int(params['pre_fac'])}_fd{int(params['fd'])}"

    def setup_run_directory(self, run_dir, params):
        """Create run directory and prepare input file"""
        run_dir = Path(run_dir)
        run_dir.mkdir(parents=True, exist_ok=True)

        # Symlink dependencies
        deps = ['initial_cells.data', 'random_propelling_directions.x']
        for dep in deps:
            src = self.nexttissue_dir / dep
            dst = run_dir / dep
            if not dst.exists() and src.exists():
                os.symlink(src.absolute(), dst)

        # Modify input file to redirect outputs to run_dir
        input_path = self.nexttissue_dir / self.input_file
        with open(input_path, 'r') as f:
            content = f.read()

        # Replace hardcoded output paths with run-specific directory
        content = content.replace('/mnt/d/NextTissUe_data/', str(run_dir) + '/')

        # Save as temporary input file
        tmp_input = run_dir / 'in.demo_hexatic.tmp'
        with open(tmp_input, 'w') as f:
            f.write(content)

        return tmp_input

    def run_simulation(self, params, run_dir, run_idx=None, total_runs=None, sweep_stats=None, simple_progress=False):
        """Execute LAMMPS with parameter overrides using MPI

        Args:
            sweep_stats: dict with 'completed', 'failed', 'skipped', 'avg_duration' for overall progress
            simple_progress: If True, use single-line progress (for terminals without ANSI support)
        """
        tmp_input = run_dir / 'in.demo_hexatic.tmp'

        # Build MPI command with parameter overrides
        cmd = [
            'mpirun', '-np', str(self.np),
            self.lammps_exe, '-in', 'in.demo_hexatic.tmp',
            # Sweep parameters (vary per run)
            '-var', 'karea', str(params['karea']),
            '-var', 'pre_fac', str(params['pre_fac']),
            '-var', 'fd', str(params['fd']),
            # Global simulation parameters
            '-var', 'main_run_steps', str(self.main_run_steps),
            '-var', 'dump_interval', str(self.dump_interval)
        ]

        self.log_progress(f"Running: mpirun -np {self.np} lmp ... (karea={params['karea']}, pre_fac={params['pre_fac']}, fd={params['fd']})")

        # Execute LAMMPS with real-time output monitoring
        try:
            process = subprocess.Popen(
                cmd,
                cwd=run_dir,
                stdout=subprocess.PIPE,
                stderr=subprocess.PIPE,
                text=True,
                bufsize=1  # Line buffered
            )

            stdout_lines = []
            stderr_lines = []
            current_step = 0
            run_start_time = datetime.now()
            last_update_time = run_start_time

            # Read output line by line for live progress
            for line in process.stdout:
                stdout_lines.append(line)

                # Parse timestep from LAMMPS thermo output
                # Look for lines like "    1000  298.456  -1234.56 ..."
                if line.strip() and not line.startswith('Step') and not line.startswith('Loop'):
                    parts = line.split()
                    if parts and parts[0].isdigit():
                        try:
                            current_step = int(parts[0])

                            # Update progress every 2 seconds to avoid spam
                            now = datetime.now()
                            if (now - last_update_time).total_seconds() >= 2.0:
                                elapsed = (now - run_start_time).total_seconds()
                                progress = current_step / self.main_run_steps

                                if progress > 0:
                                    estimated_total = elapsed / progress
                                    remaining = estimated_total - elapsed
                                    eta_run = self.format_eta(remaining)

                                    # Per-run progress bar
                                    bar = self.progress_bar(current_step, self.main_run_steps, width=30)

                                    # Line 1: Per-run progress (with compact step display)
                                    step_str = f"{self.format_number(current_step)}/{self.format_number(self.main_run_steps)}"
                                    line1 = f"  {bar} | Step {step_str} | Run ETA: {eta_run}"

                                    # Line 2: Overall sweep progress
                                    line2 = ""
                                    if sweep_stats and run_idx and total_runs:
                                        completed = sweep_stats.get('completed', 0)
                                        failed = sweep_stats.get('failed', 0)
                                        skipped = sweep_stats.get('skipped', 0)
                                        sweep_start = sweep_stats.get('sweep_start_time')

                                        # Build overall progress line
                                        runs_done = completed + failed
                                        line2 = f"  Overall: {runs_done}/{total_runs} (C:{completed} F:{failed} S:{skipped})"

                                        # Calculate overall ETA using cumulative progress (need only 0.5% progress)
                                        if sweep_start and progress > 0.005:
                                            # Total time elapsed since sweep started
                                            total_elapsed = (now - sweep_start).total_seconds()

                                            # Total progress: completed runs + current run progress
                                            total_progress = runs_done + progress

                                            # Estimate based on cumulative average speed
                                            if total_progress > 0:
                                                est_time_per_run = total_elapsed / total_progress
                                                remaining_progress = total_runs - total_progress
                                                est_remaining_seconds = remaining_progress * est_time_per_run

                                                overall_eta = self.format_eta(est_remaining_seconds)
                                                line2 += f" | ETA: {overall_eta}"

                                    # Display progress - choose method based on terminal support
                                    if simple_progress:
                                        # Fallback: Single-line display with all info (for terminals without ANSI)
                                        combined = f"{line1} | {line2}"
                                        print(f"\r{combined}", end='', flush=True)
                                    else:
                                        # Two-line display with ANSI escape codes
                                        if not hasattr(self, '_live_progress_lines'):
                                            # First time - print placeholder lines
                                            self._live_progress_lines = True
                                            print(f"{line1}")
                                            print(f"{line2}", end='', flush=True)
                                        else:
                                            # Subsequent updates - move cursor up and overwrite
                                            # Move up 1 line first, then print line1, newline, line2, and move up again
                                            print(f"\x1b[1A\r\x1b[K{line1}\n\x1b[K{line2}", end='', flush=True)

                                last_update_time = now
                        except (ValueError, IndexError):
                            pass

            # Clear the live progress lines and move to next line
            if simple_progress:
                print()  # Single line mode - just newline
            else:
                print("\n")  # Two-line mode - extra newline to clear both lines

            # Read any remaining stderr
            stderr_output = process.stderr.read()
            stderr_lines = stderr_output.split('\n') if stderr_output else []

            # Wait for process to complete
            returncode = process.wait()

            # Save stdout and stderr
            with open(run_dir / 'lammps_stdout.log', 'w') as f:
                f.writelines(stdout_lines)
            with open(run_dir / 'lammps_stderr.log', 'w') as f:
                f.writelines(stderr_lines)

            if returncode != 0:
                raise RuntimeError(f"LAMMPS exited with code {returncode}")

            return True

        except Exception as e:
            self.log_progress(f"ERROR: {str(e)}")
            return False

    def save_used_input_file(self, params, run_dir):
        """Save in.demo_hexatic.used with all ${var} substituted"""
        tmp_input = run_dir / 'in.demo_hexatic.tmp'
        used_input = run_dir / 'in.demo_hexatic.used'

        with open(tmp_input, 'r') as f:
            content = f.read()

        # Substitute all parameter variables with actual values
        substitutions = {
            '${karea}': str(params['karea']),
            '${pre_fac}': str(params['pre_fac']),
            '${fd}': str(params['fd']),
            '${main_run_steps}': str(self.main_run_steps),
            '${dump_interval}': str(self.dump_interval),
        }

        for var, val in substitutions.items():
            content = content.replace(var, val)

        with open(used_input, 'w') as f:
            f.write(content)

    def mark_complete(self, run_dir):
        """Create COMPLETE marker file"""
        marker = run_dir / 'COMPLETE'
        with open(marker, 'w') as f:
            f.write(f"Completed at {datetime.now().isoformat()}\n")

    def is_complete(self, run_dir):
        """Check if COMPLETE marker exists"""
        return (run_dir / 'COMPLETE').exists()

    def format_number(self, num):
        """Format number with k/M suffix to save space

        Examples: 15000 -> 15k, 500000 -> 500k, 1500000 -> 1.5M
        """
        if num >= 1_000_000:
            return f"{num/1_000_000:.1f}M".rstrip('0').rstrip('.')
        elif num >= 1_000:
            return f"{num/1_000:.0f}k"
        else:
            return str(num)

    def format_eta(self, seconds_remaining):
        """Format estimated time of arrival with relative day indicator

        Returns: String like "Today at 14:32" or "Tomorrow at 09:15"
        """
        finish_time = datetime.now() + timedelta(seconds=seconds_remaining)
        finish_str = finish_time.strftime('%H:%M')

        # Calculate days difference
        today = datetime.now().date()
        finish_date = finish_time.date()
        days_diff = (finish_date - today).days

        if days_diff == 0:
            return f"Today at {finish_str}"
        elif days_diff == 1:
            return f"Tomorrow at {finish_str}"
        else:
            return f"{days_diff} days later at {finish_str}"

    def progress_bar(self, current, total, width=40):
        """Generate a text progress bar

        Example: [████████████████░░░░░░░░] 62.5%
        """
        filled = int(width * current / total)
        bar = '█' * filled + '░' * (width - filled)
        percent = 100 * current / total
        return f"[{bar}] {percent:.1f}%"

    def execute_sweep(self, campaign_ids, resume=False, simple_progress=False):
        """Run all simulations sequentially

        Args:
            resume: If False (default), skip completed runs. If True, force re-run even completed runs.
            simple_progress: Use single-line progress for terminals without ANSI support
        """
        # Initialize manifest
        self.output_base.mkdir(parents=True, exist_ok=True)

        # Collect all runs
        all_runs = []
        for campaign_id in campaign_ids:
            combinations = self.generate_run_combinations(campaign_id)
            all_runs.extend(combinations)

        total_runs = len(all_runs)
        sweep_start_time = datetime.now()  # Track sweep start for cumulative ETA
        self.log_progress(f"Starting parameter sweep: {total_runs} runs across {len(campaign_ids)} campaign(s)")
        self.log_progress(f"MPI processes per run: {self.np}")
        self.log_progress(f"Main run steps: {self.main_run_steps}")
        self.log_progress(f"Dump interval: {self.dump_interval}")

        completed = 0
        skipped = 0
        failed = 0
        run_durations = []  # Track durations for ETA calculation

        for idx, params in enumerate(all_runs, 1):
            run_name = self.get_run_dir_name(params)
            run_dir = self.output_base / run_name
            campaign_id = params['campaign']

            self.log_progress(f"\n{'='*70}")
            self.log_progress(f"Run {idx}/{total_runs}: {run_name}")
            self.log_progress(f"Campaign {campaign_id} | karea={params['karea']} | pre_fac={params['pre_fac']} | fd={params['fd']}")

            # Check if already complete (always skip unless --force-rerun)
            if not resume and self.is_complete(run_dir):
                self.log_progress(f"SKIPPED: Run already complete (use --force-rerun to override)")
                self.update_manifest(params, 'completed', None, None, None, run_dir)
                skipped += 1
                completed += 1
                continue

            # Run simulation
            start_time = datetime.now()
            self.log_progress(f"Started at {start_time.strftime('%H:%M:%S')}")

            try:
                # Setup
                self.setup_run_directory(run_dir, params)

                # Prepare sweep stats for live progress display
                avg_duration = sum(run_durations) / len(run_durations) if run_durations else 0
                sweep_stats = {
                    'completed': completed,
                    'failed': failed,
                    'skipped': skipped,
                    'avg_duration': avg_duration,
                    'sweep_start_time': sweep_start_time  # For cumulative ETA calculation
                }

                # Execute with live progress
                success = self.run_simulation(params, run_dir, run_idx=idx, total_runs=total_runs,
                                              sweep_stats=sweep_stats, simple_progress=simple_progress)

                if success:
                    # Save used input file and mark complete
                    self.save_used_input_file(params, run_dir)
                    self.mark_complete(run_dir)

                    end_time = datetime.now()
                    duration = (end_time - start_time).total_seconds()
                    run_durations.append(duration)  # Track for ETA

                    self.log_progress(f"COMPLETED in {duration:.1f}s ({duration/60:.1f} min)")
                    self.update_manifest(params, 'completed', start_time, end_time, duration, run_dir)
                    completed += 1
                else:
                    end_time = datetime.now()
                    duration = (end_time - start_time).total_seconds()
                    run_durations.append(duration)  # Track even failed runs for ETA

                    self.log_progress(f"FAILED after {duration:.1f}s")
                    self.update_manifest(params, 'failed', start_time, end_time, duration, run_dir)
                    failed += 1

            except Exception as e:
                end_time = datetime.now()
                duration = (end_time - start_time).total_seconds()
                run_durations.append(duration)  # Track even errors for ETA

                self.log_progress(f"ERROR: {str(e)}")
                self.update_manifest(params, 'error', start_time, end_time, duration, run_dir)
                failed += 1

            # Progress summary with progress bar and ETA
            runs_done = completed + failed
            progress_bar = self.progress_bar(runs_done, total_runs)

            # Calculate ETA if we have duration data
            if run_durations:
                avg_duration = sum(run_durations) / len(run_durations)
                remaining_runs = total_runs - runs_done
                estimated_seconds = avg_duration * remaining_runs
                eta_str = self.format_eta(estimated_seconds)

                self.log_progress(f"{progress_bar}")
                self.log_progress(f"Progress: {runs_done}/{total_runs} | Completed: {completed} | Failed: {failed} | Skipped: {skipped}")
                self.log_progress(f"Estimated completion: {eta_str} (avg {avg_duration/60:.1f} min/run, {remaining_runs} runs left)")
            else:
                self.log_progress(f"{progress_bar}")
                self.log_progress(f"Progress: {runs_done}/{total_runs} | Completed: {completed} | Failed: {failed} | Skipped: {skipped}")

        # Final summary
        self.log_progress(f"\n{'='*70}")
        self.log_progress(f"SWEEP COMPLETE")
        self.log_progress(f"Total runs: {total_runs}")
        self.log_progress(f"Completed: {completed}")
        self.log_progress(f"Failed: {failed}")
        self.log_progress(f"Skipped (resumed): {skipped}")
        self.log_progress(f"Manifest saved to: {self.manifest_file}")

    def log_progress(self, message):
        """Log to both console and sweep.log"""
        timestamp = datetime.now().strftime('%Y-%m-%d %H:%M:%S')
        log_msg = f"[{timestamp}] {message}"

        # Print to console
        print(message)

        # Append to log file
        with open(self.log_file, 'a') as f:
            f.write(log_msg + '\n')

    def update_manifest(self, params, status, start_time, end_time, duration, run_dir):
        """Update CSV manifest with run status"""
        # Check if manifest exists
        file_exists = self.manifest_file.exists()

        # Prepare row
        row = {
            'campaign': params['campaign'],
            'karea': params['karea'],
            'pre_fac': params['pre_fac'],
            'fd': params['fd'],
            'main_run_steps': self.main_run_steps,
            'dump_interval': self.dump_interval,
            'status': status,
            'start_time': start_time.isoformat() if start_time else '',
            'end_time': end_time.isoformat() if end_time else '',
            'duration_s': f"{duration:.1f}" if duration else '',
            'output_dir': str(run_dir)
        }

        # Write to CSV
        with open(self.manifest_file, 'a', newline='') as f:
            writer = csv.DictWriter(f, fieldnames=row.keys())
            if not file_exists:
                writer.writeheader()
            writer.writerow(row)

def load_config(config_file):
    """Load configuration from YAML file

    Returns: (campaigns_dict, lammps_exe, np, main_run_steps, dump_interval)
    """
    with open(config_file, 'r') as f:
        config = yaml.safe_load(f)

    # Extract campaign definitions
    campaigns = {}
    for campaign_id, campaign_data in config.get('campaigns', {}).items():
        campaigns[campaign_id] = {
            'fixed': campaign_data.get('fixed', {}),
            'sweep': campaign_data.get('sweep', {})
        }

    # Extract execution parameters (with defaults)
    lammps_exe = config.get('lammps_exe', '/home/lost_octave/LAMMPS/build/lmp')
    np = config.get('mpi_processes', 12)

    # Extract simulation parameters (with defaults)
    sim_config = config.get('simulation', {})
    main_run_steps = sim_config.get('main_run_steps', 500000)
    dump_interval = sim_config.get('dump_interval', 5000)

    return campaigns, lammps_exe, np, main_run_steps, dump_interval

def main():
    parser = argparse.ArgumentParser(
        description='Run LAMMPS multi-parameter sweep for NextTissUe',
        formatter_class=argparse.RawDescriptionHelpFormatter,
        epilog="""
Examples:
  # Run all campaigns (32 runs) - automatically skips completed runs
  python3 run_parameter_sweep.py --all

  # Run only Campaign 1 (karea sweep)
  python3 run_parameter_sweep.py --campaign 1

  # Quick test with short run
  python3 run_parameter_sweep.py --campaign 1 --main-run-steps 10000 --dump-interval 100

  # Force re-run even if COMPLETE markers exist (rare use case)
  python3 run_parameter_sweep.py --all --force-rerun
        """
    )

    # Campaign selection
    parser.add_argument('--all', action='store_true',
                       help='Run all campaigns (1 and 2)')
    parser.add_argument('--campaign', type=int, choices=[1, 2],
                       help='Run specific campaign only')

    # Execution configuration
    parser.add_argument('--lammps-exe', default='/home/lost_octave/LAMMPS/build/lmp',
                       help='Path to LAMMPS executable')
    parser.add_argument('--np', type=int, default=12,
                       help='Number of MPI processes per simulation (default: 12)')
    parser.add_argument('--input-file', default='in.demo_hexatic',
                       help='LAMMPS input file')
    parser.add_argument('--output-base', default='/mnt/d/NextTissUe_data/run',
                       help='Base directory for outputs')

    # Simulation parameters
    parser.add_argument('--main-run-steps', type=int, default=500000,
                       help='Main measurement phase duration in steps (default: 500000, user modified)')
    parser.add_argument('--dump-interval', type=int, default=5000,
                       help='Trajectory dump frequency in steps (default: 5000, user modified)')

    # Control options
    parser.add_argument('--force-rerun', action='store_true',
                       help='Force re-run even if COMPLETE marker exists (default: skip completed runs)')
    parser.add_argument('--config', type=str,
                       help='YAML config file (if not specified, uses hardcoded campaigns)')
    parser.add_argument('--simple-progress', action='store_true',
                       help='Use single-line progress display (for terminals without ANSI support)')

    args = parser.parse_args()

    # Load configuration from YAML if provided, otherwise use hardcoded defaults
    if args.config:
        print(f"Loading configuration from {args.config}")
        config_campaigns, config_lammps_exe, config_np, config_main_run_steps, config_dump_interval = load_config(args.config)

        # Use config values as defaults, but allow command-line override
        campaigns = config_campaigns
        lammps_exe = args.lammps_exe if args.lammps_exe != parser.get_default('lammps_exe') else config_lammps_exe
        np = args.np if args.np != parser.get_default('np') else config_np
        main_run_steps = args.main_run_steps if args.main_run_steps != parser.get_default('main_run_steps') else config_main_run_steps
        dump_interval = args.dump_interval if args.dump_interval != parser.get_default('dump_interval') else config_dump_interval
    else:
        # Use hardcoded campaigns
        campaigns = {
            1: {  # Campaign 1: karea sweep
                'fixed': {'pre_fac': 10},
                'sweep': {
                    'karea': [500, 1000, 2000, 5000],
                    'fd': [200, 500, 1000, 2000]
                }
            },
            2: {  # Campaign 2: pre_fac sweep
                'fixed': {'karea': 1000},
                'sweep': {
                    'pre_fac': [5, 10, 20, 50],
                    'fd': [200, 500, 1000, 2000]
                }
            }
        }
        # Use command-line arguments
        lammps_exe = args.lammps_exe
        np = args.np
        main_run_steps = args.main_run_steps
        dump_interval = args.dump_interval

    # Determine which campaigns to run
    if args.all:
        campaign_ids = [1, 2]
    elif args.campaign:
        campaign_ids = [args.campaign]
    else:
        parser.print_help()
        print("\nError: Must specify --all or --campaign N")
        sys.exit(1)

    # Check LAMMPS executable exists
    if not Path(lammps_exe).exists():
        print(f"Error: LAMMPS executable not found: {lammps_exe}")
        sys.exit(1)

    # Create sweep instance and execute
    sweep = ParameterSweep(
        campaigns=campaigns,
        lammps_exe=lammps_exe,
        np=np,
        input_file=args.input_file,
        output_base=args.output_base,
        main_run_steps=main_run_steps,
        dump_interval=dump_interval
    )

    sweep.execute_sweep(campaign_ids, resume=args.force_rerun, simple_progress=args.simple_progress)

if __name__ == '__main__':
    main()
